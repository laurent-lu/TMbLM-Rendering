#include "tm.h"

#define TM6_LOBE_COUNT_MAX  (TM_LAYER_COUNT_MAX * 3)

MTS_NAMESPACE_BEGIN

class TM6 : public BSDF {
private:

    // Layers parameters
    int m_layer_count;
    int m_component_count;

    std::vector<ref<const Texture>> m_tex_etas;
    std::vector<ref<const Texture>> m_tex_kappas;
    std::vector<ref<const Texture>> m_tex_alphas;
    std::vector<ref<const Texture>> m_tex_depths;
    std::vector<ref<const Texture>> m_tex_sigmas_s;
    std::vector<ref<const Texture>> m_tex_sigmas_k;
    std::vector<ref<const Texture>> m_tex_gs;

    // LUTs
    lut2 m_GD;
    lut4 m_FGD;
    lut3 m_TIR;

public:

    TM6(const Properties& props): BSDF(props) {

        // Structure properties
        parseLayers(props, m_layer_count,
                    m_tex_etas, m_tex_kappas, m_tex_alphas,
                    m_tex_depths, m_tex_sigmas_s, m_tex_sigmas_k, m_tex_gs);

        m_component_count = m_layer_count * 2 - 1;

        // LUTs
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_GD    = lut2(fResolver->resolve("data/tm_GD.bin").string());
        m_FGD   = lut4(fResolver->resolve("data/tm_FGD.bin").string());
        m_TIR   = lut3(fResolver->resolve("data/tm_TIR.bin").string());
    }
    TM6(Stream* stream, InstanceManager* manager): BSDF(stream, manager) {

        configure();
    }

	void addChild(const std::string &name, ConfigurableObject *child) {

		std::string prefix = name.substr(0, name.size() - 2);
		int index = atoi(name.substr(name.size() - 1).c_str());

		if (child->getClass()->derivesFrom(MTS_CLASS(Texture2D))) {

            if (prefix == "eta") {
                SLog(EInfo, "Adding n texture for layer %d", index);
                m_tex_etas[index + 1] = static_cast<Texture2D*>(child);
			}
            else if (prefix == "kappa") {
                SLog(EInfo, "Adding k texture for layer %d", index);
                m_tex_kappas[index + 1] = static_cast<Texture2D*>(child);
			}
            else if (prefix == "alpha") {
                SLog(EInfo, "Adding alpha texture for layer %d", index);
                m_tex_alphas[index + 1] = static_cast<Texture2D*>(child);
			}

            else if (prefix == "depth") {
                SLog(EInfo, "Adding depth texture for layer %d", index);
                m_tex_depths[index + 1] = static_cast<Texture2D*>(child);
			}
            else if (prefix == "sigmas") {
                SLog(EInfo, "Adding sigma_s texture for layer %d", index);
                m_tex_sigmas_s[index + 1] = static_cast<Texture2D*>(child);
			}
            else if (prefix == "sigmak") {
                SLog(EInfo, "Adding sigma_k texture for layer %d", index);
                m_tex_sigmas_k[index + 1] = static_cast<Texture2D*>(child);
			}
            else if (prefix == "g") {
                SLog(EInfo, "Adding g texture for layer %d", index);
                m_tex_gs[index + 1] = static_cast<Texture2D*>(child);
			}
			else
				BSDF::addChild(name, child);
		}
		else {
			BSDF::addChild(name, child);
        }

        SLog(EInfo, "");
    }

    void configure() {

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

        m_usesRayDifferentials = false;

        BSDF::configure();
    }

    void    components_transfer_factors(Vector3 wi,
                                        const Spectrum* etas,
                                        const Spectrum* kappas,
                                        const Float* alphas,
                                        const Float* depths,
                                        const Spectrum* sigmas_s,
                                        const Spectrum* sigmas_k,
                                        const Float* gs,
                                        component_factors_6* ops) const {

        // Current layer attributes
        Spectrum etas_ij;

        // Layer iterations
        for (int i = 0; i < m_layer_count; ++i) {

            // IoR
            etas_ij = etas[i+1] / etas[i];

            if (kappas[i+1].isZero()) {

                // Dielectric interface
                dielectric_transfer_factors(wi, etas_ij.average(), alphas[i+1], m_FGD, ops[i*2]);

                // Propagation direction update
                wi = -ops[i*2].i.t.mean;

                // Homogeneous participating medium
                medium_transfer_factors(wi, depths[i+1], sigmas_s[i+1], sigmas_k[i+1], gs[i+1], ops[i*2+1]);
            }
            else {

                // Conductor substrate
                conductor_transfer_factors(wi, etas_ij, kappas[i+1] / etas[i], alphas[i+1], m_FGD, ops[i*2]);
            }
        }
    }
    int     outgoing_lobes(const Vector3& wi,
                           const Spectrum* etas,
                           const Spectrum* kappas,
                           const Float* alphas,
                           const Float* depths,
                           const Spectrum* sigmas_s,
                           const Spectrum* sigmas_k,
                           const Float* gs,
                           hg* lobes,
                           Vector3* lobes_wi) const {

        // Mirror direction
        const Vector3 wi_r = reflectZ(wi);

        // Current lobe index
        int l = 0;

        // Components transfer factors
        component_factors_6 ops[TM6_LOBE_COUNT_MAX];

        // Current component attributes
        Float eta_ij, eta_ji;
        Spectrum e_r_i;

        // Cumulated components attributes
        Spectrum e_r_0h(0.f), e_r_0i(0.f);
        Spectrum e_r_f_0h(0.f), e_r_f_0i(0.f);
        Spectrum e_r_b_0h(0.f), e_r_b_0i(0.f);

        Spectrum g_r_0h = Spectrum(0.f), g_r_0i = Spectrum(0.f);
        Spectrum g_r_f_0h = Spectrum(0.f), g_r_f_0i = Spectrum(0.f);
        Spectrum g_r_b_0h = Spectrum(0.f), g_r_b_0i = Spectrum(0.f);

        // First-order terms
        Float g_T_0i = 1.f, g_T_0j = 1.f, g_T_0j_R = 1.f, g_T_0j_RT = 1.f;
        Float g_s_T_0i = 1.f, g_s_T_0j = 1.f, g_s_T_0j_R = 1.f, g_s_T_0j_RT = 1.f;

        // TIR
        Spectrum tir_norm;

        // Transfer matrices
        unsigned char req = ENone;

        transfer_matrix_6s etm_0i = transfer_matrix_6s::identity();
        transfer_matrix_6s gtm_0i = transfer_matrix_6s::identity();

        ////////////////////////////////////////////////////////////////////////////////
        //
        // For readibility purpose, we first extract the transfer factors of each
        // component before the actual transfer matrix calculus loop
        //
        ////////////////////////////////////////////////////////////////////////////////

        components_transfer_factors(wi, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs, ops);

        ////////////////////////////////////////////////////////////////////////////////
        //
        // Components iterations
        //
        //  1) Compute first-order dielectric transmission factors and account for TIR
        //  2) Compute transfer matrices
        //  3) Compute outgoing lobe statistics based upon the transfer matrices
        //
        ////////////////////////////////////////////////////////////////////////////////

        for (int c = 0; c < m_component_count; ++c) {

            if (ops[c].type == ENoComponent) {

                lobes[l].norm   = Spectrum(0.f);
                lobes[l].g      = 0.f;
                lobes_wi[l++]   = wi;
                continue;
            }
            else if (ops[c].type == EDielectricInterface) {

                // Layer index
                const int i = c >> 1;

                // Scalar IoR
                eta_ij = (etas[i+1] / etas[i]).average();
                eta_ji = 1.f / eta_ij;

                // First-order terms
                g_T_0j      = hg_refract(g_T_0i, eta_ji) * ops[c].i.t.g;
                g_T_0j_R    = g_T_0j * ops[c+2].i.r.g;
                g_T_0j_RT   = hg_refract(g_T_0j_R, eta_ij) * ops[c].i.tp.g;

                g_s_T_0j    = hg_refract(g_s_T_0i, eta_ji) * ops[c].i.t_ff.g;
                g_s_T_0j_R  = g_s_T_0j * ops[c+1].m.t_ff.g * ops[c+2].i.r_ff.g * ops[c+1].m.t_ff.g;
                g_s_T_0j_RT = hg_refract(g_s_T_0j_R, eta_ij) * ops[c].i.tp_ff.g;

                // Downward transmission
                ops[c].i.t.g       = g_T_0i != 0.f ? g_T_0j / g_T_0i : 0.f;
                ops[c].i.t_ff.g    = g_s_T_0i != 0.f ? g_s_T_0j / g_s_T_0i : 0.f;

                // Upward transmission
                ops[c].i.tp.g      = g_T_0j_R != 0.f ? g_T_0j_RT / g_T_0j_R : 0.f;
                ops[c].i.tp_ff.g   = g_s_T_0j_R != 0.f ? g_s_T_0j_RT / g_s_T_0j_R : 0.f;

                if (eta_ij < 1.f) {

                    // Primary downward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[c].i.r.mean.z), hg_to_ggx(g_T_0i), eta_ij) * ops[c].i.t.norm;

                    ops[c].i.r.norm += tir_norm;
                    ops[c].i.t.norm -= tir_norm;

                    // Secondary downward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[c].i.r_ff.mean.z), hg_to_ggx(g_s_T_0i), eta_ij) * ops[c].i.t_ff.norm;

                    ops[c].i.r_ff.norm += tir_norm;
                    ops[c].i.t_ff.norm -= tir_norm;
                }
                else {

                    // Primary upward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[c].i.t.mean.z), hg_to_ggx(g_T_0j_R), eta_ji) * ops[c].i.tp.norm;

                    ops[c].i.rp.norm += tir_norm;
                    ops[c].i.tp.norm -= tir_norm;

                    // Secondary upward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[c].i.t_ff.mean.z), hg_to_ggx(g_s_T_0j_R), eta_ji) * ops[c].i.tp_ff.norm;

                    ops[c].i.rp_ff.norm += tir_norm;
                    ops[c].i.tp_ff.norm -= tir_norm;
                }

                ////////////////////////////////////////////////////////////////
                // Transfer matrices evaluation
                ////////////////////////////////////////////////////////////////

                req |= EPrimaryReflectance;
                req |= req & ESecondaryBackwardReflectance ? ESecondaryForwardReflectance : ENone;

                etm_0i *= energy_matrix(ops[c]);
                gtm_0i *= asymmetry_matrix(ops[c]);

                etm_0i.factors(e_r_0i, e_r_f_0i, e_r_b_0i, req);
                gtm_0i.factors(g_r_0i, g_r_f_0i, g_r_b_0i, req);
            }
            else if (ops[c].type == EHomogeneousMedium) {

                ////////////////////////////////////////////////////////////////
                // Transfer matrices evaluation
                ////////////////////////////////////////////////////////////////

                req |= ESecondaryBackwardReflectance;
                req |= req & EPrimaryReflectance ? ESecondaryForwardReflectance : ENone;

                etm_0i *= energy_matrix(ops[c]);
                gtm_0i *= asymmetry_matrix(ops[c]);

                etm_0i.factors(e_r_0i, e_r_f_0i, e_r_b_0i, req);
                gtm_0i.factors(g_r_0i, g_r_f_0i, g_r_b_0i, req);
            }
            else {

                ////////////////////////////////////////////////////////////////
                // Transfer matrices evaluation
                ////////////////////////////////////////////////////////////////

                req |= EPrimaryReflectance;
                req |= req & ESecondaryBackwardReflectance ? ESecondaryForwardReflectance : ENone;

                etm_0i.factors(e_r_0i, e_r_f_0i, e_r_b_0i, ops[c].i.r.norm, req);
                gtm_0i.factors(g_r_0i, g_r_f_0i, g_r_b_0i, ops[c].i.r.norm * ops[c].i.r.g, req);
            }

            ////////////////////////////////////////////////////////////////////
            // Primary reflection lobe
            ////////////////////////////////////////////////////////////////////

            if (req & EPrimaryReflectance) {

                // Energy
                e_r_i = e_r_0i - e_r_0h;

                // Lobe attributes
                lobes[l].norm   = max(e_r_i, 0.f);
                lobes[l].g      = safe_div(g_r_0i - g_r_0h, e_r_i).average();
                lobes_wi[l++] = wi;
            }

            ////////////////////////////////////////////////////////////////////
            // Secondary forward reflection lobe
            ////////////////////////////////////////////////////////////////////

            if (req & ESecondaryForwardReflectance) {

                // Energy
                e_r_i = e_r_f_0i - e_r_f_0h;

                // Lobe attributes
                lobes[l].norm   = max(e_r_i, 0.f);
                lobes[l].g      = safe_div(g_r_f_0i - g_r_f_0h, e_r_i).average();
                lobes_wi[l++] = wi;
            }

            ////////////////////////////////////////////////////////////////////
            // Secondary backward reflection lobe
            ////////////////////////////////////////////////////////////////////

            if (req & ESecondaryBackwardReflectance) {

                // Energy
                e_r_i = e_r_b_0i - e_r_b_0h;

                // Lobe attributes
                lobes[l].norm   = max(e_r_i, 0.f);
                lobes[l].g      = safe_div(g_r_b_0i - g_r_b_0h, e_r_i).average();
                lobes_wi[l++] = wi_r;
            }

            // Top layers attributes update
            e_r_0h = e_r_0i;
            g_r_0h = g_r_0i;

            e_r_f_0h = e_r_f_0i;
            g_r_f_0h = g_r_f_0i;

            e_r_b_0h = e_r_b_0i;
            g_r_b_0h = g_r_b_0i;

            // First-order terms update
            g_T_0i      = g_T_0j;
            g_s_T_0i    = g_s_T_0j;
        }

        return l;
    }

    inline void evalMaps(const BSDFSamplingRecord &bRec,
                         Spectrum* etas, Spectrum* kappas, Float* alphas,
                         Float* depths, Spectrum* sigmas_s, Spectrum* sigmas_k, Float* gs) const {

        // External medium properties
        etas[0] = m_tex_etas[0]->eval(bRec.its);

        // Layers properties
        for(int i = 1; i <= m_layer_count; ++i) {

            etas[i]   = m_tex_etas[i]->eval(bRec.its);
            kappas[i] = m_tex_kappas[i]->eval(bRec.its);
            alphas[i] = m_tex_alphas[i]->eval(bRec.its).average();

            depths[i] = m_tex_depths[i]->eval(bRec.its).average();
            sigmas_s[i] = m_tex_sigmas_s[i]->eval(bRec.its);
            sigmas_k[i] = m_tex_sigmas_k[i]->eval(bRec.its); sigmas_k[i] = max(sigmas_k[i], 1e-7f);
            gs[i] = m_tex_gs[i]->eval(bRec.its).average();
        }
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        Spectrum    kappas[TM_LAYER_COUNT_MAX];
        Float       alphas[TM_LAYER_COUNT_MAX];
        Float       depths[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_s[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_k[TM_LAYER_COUNT_MAX];
        Float       gs[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs);

        return this->eval(bRec, measure, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs);
    }
    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure,
                  const Spectrum* etas, const Spectrum* kappas, const Float* alphas,
                  const Float* depths, const Spectrum* sigmas_s, const Spectrum* sigmas_k, const Float* gs) const {

        if (measure != ESolidAngle || bRec.wi.z <= 0 || bRec.wo.z <= 0)
            return Spectrum(0.0f);

        // Half-vector
        Vector3 H = normalize(bRec.wi + bRec.wo);

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        hg      lobes[TM6_LOBE_COUNT_MAX];
        Vector3 lobes_wi[TM6_LOBE_COUNT_MAX];

        const int lobe_count = outgoing_lobes(bRec.wi,
                                              etas, kappas, alphas,
                                              depths, sigmas_s, sigmas_k, gs,
                                              lobes, lobes_wi);

        ////////////////////////////////////////////////////////////////////////
        // Throughput
        ////////////////////////////////////////////////////////////////////////

        Spectrum throughput = Spectrum(0.f);

        // [Bati2019] Top reflection correction
        if (!lobes[0].norm.isZero()) {

            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, alphas[1], true);

            const Float G2  = ndf.G(bRec.wi, bRec.wo, H);
            const Float D   = ndf.eval(H);

            Spectrum F;
            const Spectrum etas_01 = etas[1] / etas[0];
            if (kappas[1].isZero())
                F = Spectrum(fresnelDielectricExt(dot(bRec.wi, H), etas_01.average()));
            else
                F = fresnelConductorExact(dot(bRec.wi, H), etas_01, kappas[1] / etas[0]);

            throughput += F * G2 * D / (4.f * bRec.wi.z);
        }

        // Eval and sum internal lobes contributions
        for (int i = 1; i < lobe_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            // Incident direction
            const Vector3 wi = lobes_wi[i];

            // Half-vector
            H = normalize(wi + bRec.wo);

            // GGX equivalent normal distribution
            const Float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const Float G2  = ndf.G(wi, bRec.wo, H);
            const Float D   = ndf.eval(H);

            // Single-scattering albedo normalization
            const Float essi = m_GD.range_get_interpolate(bRec.wi.z, a);

            // Ideal BRDF eval
            const Float f = G2 * D / (4.f * wi.z * essi);

            throughput += lobes[i].norm * f;
        }

        Assert(!throughput.isNaN());

        return throughput;
    }

    Float pdf(const BSDFSamplingRecord& bRec, EMeasure measure) const {

        if (measure != ESolidAngle || bRec.wi.z <= 0 || bRec.wo.z <= 0)
            return 0.f;

        // Half-vector
        Vector3 H = normalize(bRec.wi + bRec.wo);

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        Spectrum    kappas[TM_LAYER_COUNT_MAX];
        Float       alphas[TM_LAYER_COUNT_MAX];
        Float       depths[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_s[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_k[TM_LAYER_COUNT_MAX];
        Float       gs[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs);

        hg      lobes[TM6_LOBE_COUNT_MAX];
        Vector3 lobes_wi[TM6_LOBE_COUNT_MAX];

        const int lobe_count = outgoing_lobes(bRec.wi,
                                              etas, kappas, alphas,
                                              depths, sigmas_s, sigmas_k, gs,
                                              lobes, lobes_wi);

        ////////////////////////////////////////////////////////////////////////
        // PDF
        ////////////////////////////////////////////////////////////////////////

        Float w_sum = 0.f;
        Float wpdf_sum = 0.f;
        for (int i = 0; i < lobe_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            // Incident direction
            const Vector3 wi = lobes_wi[i];

            // Half-vector
            H = normalize(wi + bRec.wo);

            // GGX equivalent normal distribution
            const Float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const Float G1  = ndf.smithG1(wi, H);
            const Float D   = ndf.eval(H);

            const Float norm = lobes[i].norm.average();

            w_sum       += norm;
            wpdf_sum    += norm * G1 * D / (4.0f * wi.z);
        }

        Assert(!(std::isnan(wpdf_sum) || std::isinf(wpdf_sum)));

        return w_sum > 0.f ? wpdf_sum / w_sum : 0.f;
    }

    Spectrum sample(BSDFSamplingRecord& bRec, const Point2& sample) const {

        Float dummy;
        return this->sample(bRec, dummy, sample);
    }
    Spectrum sample(BSDFSamplingRecord& bRec, Float& pdf, const Point2& sample) const {

        if (bRec.wi.z < 0)
            return Spectrum(0.f);

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        Spectrum    kappas[TM_LAYER_COUNT_MAX];
        Float       alphas[TM_LAYER_COUNT_MAX];
        Float       depths[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_s[TM_LAYER_COUNT_MAX];
        Spectrum    sigmas_k[TM_LAYER_COUNT_MAX];
        Float       gs[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs);

        hg      lobes[TM6_LOBE_COUNT_MAX];
        Vector3 lobes_wi[TM6_LOBE_COUNT_MAX];

        const int lobe_count = outgoing_lobes(bRec.wi,
                                              etas, kappas, alphas,
                                              depths, sigmas_s, sigmas_k, gs,
                                              lobes, lobes_wi);

        //////////////////////////////////////////////////////////////////////
        // Random lobe selection
        //////////////////////////////////////////////////////////////////////

        // Lobes probabilities
        Float w[TM6_LOBE_COUNT_MAX];
        Float w_sum = 0.f;
        for(int i = 0; i < lobe_count; ++i) {

            w[i] = lobes[i].norm.average();
            w_sum += w[i];
        }

        // Lobe selection
        Float sel_w = bRec.sampler->next1D() * w_sum - w[0];
        int   sel_i = 0;
        for(sel_i = 0; sel_w > 0.f && sel_i < lobe_count; ++sel_i)
            sel_w -= w[sel_i + 1];

        //////////////////////////////////////////////////////////////////////
        // Sampling
        //////////////////////////////////////////////////////////////////////

        // Selected lobe equivalent GGX normal distribution
        const Float sel_a = hg_to_ggx(lobes[sel_i].g);
        const MicrofacetDistribution sel_ndf(MicrofacetDistribution::EGGX, sel_a, true);

        // Selected lobe incident direction
        Vector3 sel_wi = lobes_wi[sel_i];

        // Selected lobe normal sampling
        Vector3 H = sel_ndf.sample(sel_wi, sample, pdf);

        // Sampled direction
        bRec.wo = reflect(sel_wi, H);
        bRec.eta = 1.f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        if(bRec.wo.z <= 0.f || pdf <= 0.f)
            return Spectrum(0.0f);

        //////////////////////////////////////////////////////////////////////
        // PDF
        //////////////////////////////////////////////////////////////////////

        pdf = 0.f;
        for(int i = 0; i < lobe_count; ++i) {

            if (w[i] == 0.f)
                continue;

            // Incident direction
            Vector3 wi = lobes_wi[i];

            // Half-vector
            H = normalize(wi + bRec.wo);

            // GGX equivalent normal distribution
            const Float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const Float G1 = ndf.smithG1(wi, H);
            const Float D = ndf.eval(H);

            pdf += w[i] * G1 * D / (4.0f * wi.z);
        }
        pdf /= w_sum;

        //////////////////////////////////////////////////////////////////////
        // Throughput
        //////////////////////////////////////////////////////////////////////

        Spectrum throughput = this->eval(bRec, ESolidAngle, etas, kappas, alphas, depths, sigmas_s, sigmas_k, gs);

        Assert(!throughput.isNaN());

        return pdf > 0.f ? throughput / pdf : Spectrum(0.0f);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(TM6, false, BSDF)
MTS_EXPORT_PLUGIN(TM6, "TM6");
MTS_NAMESPACE_END

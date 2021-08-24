#include "tm.h"

MTS_NAMESPACE_BEGIN

class TM2 : public BSDF {
private:

    // Layers parameters
    int m_layer_count;

    std::vector<ref<const Texture>> m_tex_etas;
    std::vector<ref<const Texture>> m_tex_kappas;
    std::vector<ref<const Texture>> m_tex_alphas;

    // LUTs
    lut4 m_FGD;
    lut3 m_TIR;

public:

    TM2(const Properties& props): BSDF(props) {

        // Structure properties
        parseLayers(props, m_layer_count, m_tex_etas, m_tex_kappas, m_tex_alphas, true);

        // LUTs
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_FGD = lut4(fResolver->resolve("data/tm_FGD.bin").string());
        m_TIR = lut3(fResolver->resolve("data/tm_TIR.bin").string());
    }
    TM2(Stream* stream, InstanceManager* manager): BSDF(stream, manager) {

        configure();
    }

    void configure() {

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | ESpatiallyVarying);

        m_usesRayDifferentials = false;

        BSDF::configure();
    }

    void components_transfer_factors(Vector3 wi,
                                     const Spectrum* etas,
                                     const Spectrum* kappas,
                                     const Float*    alphas,
                                     component_factors_2* ops) const {

        // Current layer attributes
        Spectrum etas_ij;

        // Layer iterations
        for (int i = 0; i < m_layer_count; ++i) {

            // IoR
            etas_ij = etas[i+1] / etas[i];

            if (kappas[i+1].isZero()) {

                // Dielectric interface
                dielectric_transfer_factors(wi, etas_ij.average(), alphas[i+1], m_FGD, ops[i]);

                // Propagation direction update
                wi = -ops[i].t.mean;
            }
            else {

                // Conductor substrate
                conductor_transfer_factors(wi, etas_ij, kappas[i+1] / etas[i], alphas[i+1], m_FGD, ops[i]);
            }
        }
    }
    void outgoing_lobes(const Vector3& wi,
                        const Spectrum* etas,
                        const Spectrum* kappas,
                        const Float* alphas,
                        hg* lobes) const {

        // Components transfer factors
        component_factors_2 ops[TM_LAYER_COUNT_MAX];

        // Current component attributes
        Float eta_ij, eta_ji;
        Spectrum e_r_i;
        Float e_r_i_avg;

        // Cumulated components attributes
        Spectrum e_r_0h(0.f), e_r_0i(0.f);
        Float g_r_0h = 0.f, g_r_0i;

        // First-order terms
        Float g_T_0i = 1.f, g_T_0j = 1.f, g_T_0j_R = 1.f, g_T_0j_RT = 1.f;

        // TIR
        Spectrum tir_norm;

        // Transfer matrices
        transfer_matrix_2s etm_0i = transfer_matrix_2s::identity();
        transfer_matrix_21 gtm_0i = transfer_matrix_21::identity();

        ////////////////////////////////////////////////////////////////////////////////
        //
        // For readibility purpose, we first extract the transfer factors of each
        // component before the actual transfer matrix calculus loop
        //
        ////////////////////////////////////////////////////////////////////////////////

        components_transfer_factors(wi, etas, kappas, alphas, ops);

        ////////////////////////////////////////////////////////////////////////////////
        //
        // Components iterations
        //
        //  1) Compute first-order dielectric transmission factors and account for TIR
        //  2) Compute transfer matrices
        //  3) Compute outgoing lobe statistics based upon the transfer matrices
        //
        ////////////////////////////////////////////////////////////////////////////////

        for (int i = 0; i < m_layer_count; ++i) {

            if (ops[i].type == EDielectricInterface) {

                // Scalar IoR
                eta_ij = (etas[i+1] / etas[i]).average();
                eta_ji = 1.f / eta_ij;

                // First-order terms
                g_T_0j      = hg_refract(g_T_0i, eta_ji) * ops[i].t.g;
                g_T_0j_R    = g_T_0j * ops[i+1].r.g;
                g_T_0j_RT   = hg_refract(g_T_0j_R, eta_ij) * ops[i].tp.g;

                // Downward transmission
                ops[i].t.g = g_T_0i != 0.f ? g_T_0j / g_T_0i : 0.f;

                // Upward transmission
                ops[i].tp.g = g_T_0j_R != 0.f ? g_T_0j_RT / g_T_0j_R : 0.f;

                if (eta_ij < 1.f) {

                    // Downward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[i].r.mean.z), hg_to_ggx(g_T_0i), eta_ij) * ops[i].t.norm;

                    ops[i].r.norm += tir_norm;
                    ops[i].t.norm -= tir_norm;
                }
                else {

                    // Upward TIR
                    tir_norm = m_TIR.range_get_interpolate(std::abs(ops[i].t.mean.z), hg_to_ggx(g_T_0j_R), eta_ji) * ops[i].tp.norm;

                    ops[i].rp.norm += tir_norm;
                    ops[i].tp.norm -= tir_norm;
                }

                ////////////////////////////////////////////////////////////////
                // Transfer matrices evaluation
                ////////////////////////////////////////////////////////////////
                
                etm_0i *= energy_matrix(ops[i]);
                gtm_0i *= asymmetry_matrix(ops[i]);

                e_r_0i = etm_0i.r();
                g_r_0i = gtm_0i.r();
            }
            else {

                ////////////////////////////////////////////////////////////////
                // Transfer matrices evaluation
                ////////////////////////////////////////////////////////////////

                e_r_0i = etm_0i.r(ops[i].r.norm);
                g_r_0i = gtm_0i.r(ops[i].r.norm.average() * ops[i].r.g);
            }

            ////////////////////////////////////////////////////////////////////
            // Reflected lobe
            ////////////////////////////////////////////////////////////////////

            // Energy
            e_r_i       = e_r_0i - e_r_0h;
            e_r_i_avg   = e_r_i.average();

            // Lobe attributes
            lobes[i].norm   = e_r_i;
            lobes[i].g      = e_r_i_avg > 0.f ? std::min((g_r_0i - g_r_0h) / e_r_i_avg, 1.f) : 0.f;

            // Top layers attributes update
            e_r_0h = e_r_0i;
            g_r_0h = g_r_0i;

            // First-order terms update
            g_T_0i = g_T_0j;
        }
    }

    inline void evalMaps(const BSDFSamplingRecord &bRec,
                         Spectrum* etas, Spectrum* kappas, Float* alphas) const {

        // External medium properties
        etas[0] = m_tex_etas[0]->eval(bRec.its);

        // Layers properties
        for(int i = 1; i <= m_layer_count; ++i) {

            etas[i]   = m_tex_etas[i]->eval(bRec.its);
            kappas[i] = m_tex_kappas[i]->eval(bRec.its);
            alphas[i] = m_tex_alphas[i]->eval(bRec.its).average();
        }
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        Spectrum    kappas[TM_LAYER_COUNT_MAX];
        Float       alphas[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, kappas, alphas);

        return this->eval(bRec, measure, etas, kappas, alphas);
    }
    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure,
                  const Spectrum* etas, const Spectrum* kappas, const Float* alphas) const {

        if (measure != ESolidAngle || bRec.wi.z <= 0 || bRec.wo.z <= 0)
            return Spectrum(0.0f);

        // Half-vector
        const Vector3 H = normalize(bRec.wi + bRec.wo);

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        hg lobes[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, kappas, alphas, lobes);

        ////////////////////////////////////////////////////////////////////////
        // Throughput
        ////////////////////////////////////////////////////////////////////////

        // [Bati2019] Top reflection correction
        const MicrofacetDistribution ndf_0(MicrofacetDistribution::EGGX, alphas[1], true);

        const Float G2_0    = ndf_0.G(bRec.wi, bRec.wo, H);
        const Float D_0     = ndf_0.eval(H);

        Spectrum F_0;
        const Spectrum etas_01 = etas[1] / etas[0];
        if (kappas[1].isZero())
            F_0 = Spectrum(fresnelDielectricExt(dot(bRec.wi, H), etas_01.average()));
        else
            F_0 = fresnelConductorExact(dot(bRec.wi, H), etas_01, kappas[1] / etas[0]);

        Spectrum throughput = F_0 * G2_0 * D_0 / (4.f * bRec.wi.z);

        // Eval and sum internal lobes contributions
        for (int i = 1; i < m_layer_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            // GGX equivalent normal distribution
            const Float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const Float G2  = ndf.G(bRec.wi, bRec.wo, H);
            const Float D   = ndf.eval(H);

            // Ideal BRDF eval
            const Float f = G2 * D / (4.f * bRec.wi.z);

            throughput += lobes[i].norm * f;
        }

        Assert(!throughput.isNaN());

        return throughput;
    }

    Float pdf(const BSDFSamplingRecord& bRec, EMeasure measure) const {

        if (measure != ESolidAngle || bRec.wi.z <= 0 || bRec.wo.z <= 0)
            return 0.f;

        // Half-vector
        const Vector3 H = normalize(bRec.wi + bRec.wo);

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        Spectrum    kappas[TM_LAYER_COUNT_MAX];
        Float       alphas[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, kappas, alphas);

        hg lobes[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, kappas, alphas, lobes);

        ////////////////////////////////////////////////////////////////////////
        // PDF
        ////////////////////////////////////////////////////////////////////////

        Float w_sum = 0.f;
        Float wpdf_sum = 0.f;
        for (int i = 0; i < m_layer_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            // GGX equivalent normal distribution
            const Float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const Float G1  = ndf.smithG1(bRec.wi, H);
            const Float D   = ndf.eval(H);

            const Float norm = lobes[i].norm.average();

            w_sum       += norm;
            wpdf_sum    += norm * G1 * D / (4.0f * bRec.wi.z);
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

        evalMaps(bRec, etas, kappas, alphas);

        hg lobes[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, kappas, alphas, lobes);

        //////////////////////////////////////////////////////////////////////
        // Random lobe selection
        //////////////////////////////////////////////////////////////////////

        // Lobes probabilities
        Float w[TM_LAYER_COUNT_MAX];
        Float w_sum = 0.f;
        for(int i = 0; i < m_layer_count; ++i) {

            w[i] = lobes[i].norm.average();
            w_sum += w[i];
        }

        // Lobe selection
        Float sel_w = bRec.sampler->next1D() * w_sum - w[0];
        int   sel_i = 0;
        for(sel_i = 0; sel_w > 0.f && sel_i < m_layer_count; ++sel_i)
            sel_w -= w[sel_i + 1];

        //////////////////////////////////////////////////////////////////////
        // Sampling
        //////////////////////////////////////////////////////////////////////

        // Selected lobe equivalent GGX normal distribution
        const Float sel_a = hg_to_ggx(lobes[sel_i].g);
        const MicrofacetDistribution sel_ndf(MicrofacetDistribution::EGGX, sel_a, true);

        // Selected lobe normal sampling
        const Normal H = sel_ndf.sample(bRec.wi, sample, pdf);

        // Sampled direction
        bRec.wo = reflect(bRec.wi, H);
        bRec.eta = 1.f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        if(bRec.wo.z <= 0.f || pdf <= 0.f)
            return Spectrum(0.0f);

        //////////////////////////////////////////////////////////////////////
        // PDF
        //////////////////////////////////////////////////////////////////////

        pdf = 0.f;
        for(int i = 0; i < m_layer_count; ++i) {

            if (w[i] > 0.f) {

                // GGX equivalent normal distribution
                const Float a = hg_to_ggx(lobes[i].g);
                const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

                const Float G1  = ndf.smithG1(bRec.wi, H);
                const Float D   = ndf.eval(H);

                pdf += w[i] * G1 * D / (4.0f * bRec.wi.z);
            }
        }
        pdf /= w_sum;

        //////////////////////////////////////////////////////////////////////
        // Throughput
        //////////////////////////////////////////////////////////////////////

        Spectrum throughput = this->eval(bRec, ESolidAngle, etas, kappas, alphas);

        Assert(!throughput.isNaN());

        return pdf > 0.f ? throughput / pdf : Spectrum(0.0f);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(TM2, false, BSDF)
MTS_EXPORT_PLUGIN(TM2, "TM2");
MTS_NAMESPACE_END

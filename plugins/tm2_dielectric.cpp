#include "tm.h"

MTS_NAMESPACE_BEGIN

class TM2_dielectric : public BSDF {
private:

    // Layers parameters
    int m_layer_count;
    int m_lobe_count;

    std::vector<ref<const Texture>> m_tex_etas;
    std::vector<ref<const Texture>> m_tex_kappas;
    std::vector<ref<const Texture>> m_tex_alphas;

    // LUTs
    lut4 m_FGD;
    lut3 m_TIR;

public:

    TM2_dielectric(const Properties& props): BSDF(props) {

        // Structure properties
        parseLayers(props, m_layer_count, m_tex_etas, m_tex_kappas, m_tex_alphas, false);

        m_lobe_count = m_layer_count + 1;

        // LUTs
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_FGD = lut4(fResolver->resolve("data/tm_FGD.bin").string());
        m_TIR = lut3(fResolver->resolve("data/tm_TIR.bin").string());
    }
    TM2_dielectric(Stream* stream, InstanceManager* manager): BSDF(stream, manager) {

        configure();
    }

    void configure() {

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | EBackSide | ESpatiallyVarying);
        m_components.push_back(EGlossyTransmission | EFrontSide | EBackSide | ESpatiallyVarying | ENonSymmetric);

        m_usesRayDifferentials = false;

        BSDF::configure();
    }

    void components_transfer_factors(Vector3 wi,
                                     const Spectrum* etas,
                                     const float*    alphas,
                                     component_factors_2* ops) const {

        // Current layer attributes
        Spectrum etas_ij;

        // Layer iterations
        for (int i = 0; i < m_layer_count; ++i) {

            // Dielectric interface
            dielectric_transfer_factors(wi, (etas[i+1] / etas[i]).average(), alphas[i+1], m_FGD, ops[i]);

            // Propagation direction update
            wi = -ops[i].t.mean;
        }
    }
    void outgoing_lobes(const Vector3& wi,
                        const Spectrum* etas,
                        const float* alphas,
                        hg* lobes,
                        Vector3* lobes_wi) const {

        // Components transfer factors
        component_factors_2 ops[TM_LAYER_COUNT_MAX];

        // Current interface attributes
        float eta_ij, eta_ji;
        Spectrum e_r_i;
        float e_r_i_avg;

        // Cumulated interfaces attributes
        Spectrum e_r_0h(0.f), e_r_0i, e_t_0i;
        float e_t_0i_avg;
        float g_r_0h = 0.f, g_r_0i, g_t_0i;

        // First-order terms
        float g_T_0i = 1.f, g_T_0j, g_T_0j_R, g_T_0j_RT;

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

        components_transfer_factors(wi, etas, alphas, ops);

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

            // IoR
            eta_ij = (etas[i+1] / etas[i]).average();
            eta_ji = 1.f / eta_ij;

            ////////////////////////////////////////////////////////////////////
            // Transmission correction
            ////////////////////////////////////////////////////////////////////

            // First-order terms
            g_T_0j      = hg_refract(g_T_0i, eta_ji) * ops[i].t.g;
            g_T_0j_R    = i < (m_layer_count - 1) ? g_T_0j * ops[i+1].r.g : 1.f;
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

            ////////////////////////////////////////////////////////////////////
            // Transfer matrices evaluation
            ////////////////////////////////////////////////////////////////////

            etm_0i *= energy_matrix(ops[i]);
            gtm_0i *= asymmetry_matrix(ops[i]);

            ////////////////////////////////////////////////////////////////////
            // Reflected lobe
            ////////////////////////////////////////////////////////////////////

            // Norm
            e_r_0i      = etm_0i.r();
            e_r_i       = e_r_0i - e_r_0h;
            e_r_i_avg   = e_r_i.average();

            // Asymmetry parameter
            g_r_0i = gtm_0i.r();

            // Reflected lobe
            lobes[i].norm   = e_r_i;
            lobes[i].g      = e_r_i_avg > 0.f ? std::min((g_r_0i - g_r_0h) / e_r_i_avg, 1.f) : 0.f;
            lobes_wi[i] = wi;

            // Top layers attributes update
            e_r_0h = e_r_0i;
            g_r_0h = g_r_0i;

            // First-order terms update
            g_T_0i = g_T_0j;
        }

        ////////////////////////////////////////////////////////////////////////
        // Transmitted lobe
        ////////////////////////////////////////////////////////////////////////

        // Norm
        e_t_0i      = etm_0i.t();
        e_t_0i_avg  = e_t_0i.average();

        // Asymmetry parameter
        g_t_0i = gtm_0i.t();

        // Transmitted lobe
        lobes[m_layer_count].norm   = e_t_0i;
        lobes[m_layer_count].g      = e_t_0i_avg > 0.f ? std::min(g_t_0i / e_t_0i_avg, 1.f) : 0.f;
        lobes_wi[m_layer_count] = reflectZ(ops[m_layer_count - 1].t.mean);
    }

    inline void evalMaps(const BSDFSamplingRecord &bRec, Spectrum* etas, float* alphas) const {

        const bool rev = bRec.wi.z < 0.f;

        // External medium properties
        etas[0] = m_tex_etas[rev ? m_layer_count : 0]->eval(bRec.its);

        // Layers properties
        for(int i = 1; i <= m_layer_count; ++i) {

            etas[i]   = m_tex_etas[rev ? m_layer_count - i : i]->eval(bRec.its);
            alphas[i] = m_tex_alphas[rev ? m_layer_count - i + 1 : i]->eval(bRec.its).average();
        }
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        float       alphas[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, alphas);

        return this->eval(bRec, measure, etas, alphas);
    }
    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure,
                  const Spectrum* etas, const float* alphas) const {

        if (measure != ESolidAngle || Frame::cosTheta(bRec.wi) == 0)
            return Spectrum(0.0f);

        // Interaction type
        bool reflect = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0;

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        hg      lobes[TM_LAYER_COUNT_MAX];
        Vector3 lobes_wi[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, alphas, lobes, lobes_wi);

        ////////////////////////////////////////////////////////////////////////
        // Throughput
        ////////////////////////////////////////////////////////////////////////

        // Radiance transmission scaling factor
        const float factor = reflect ? 1.f : (etas[0] / etas[m_layer_count]).average();

        // Eval and sum lobe contributions according to the interaction type
        Spectrum throughput = Spectrum(0.f);
        for (int i = 0; i < m_lobe_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            if (lobes_wi[i].z * bRec.wo.z <= 0.f)
                continue;

            // Eval directions (upper hemisphere)
            Vector3 wi = lobes_wi[i];
            Vector3 wo = bRec.wo;

            wi.z = std::abs(wi.z);
            wo.z = std::abs(wo.z);

            // Half-vector (upper hemisphere)
            const Vector3 H = normalize(wi + wo);

            // GGX equivalent normal distribution
            const float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const float G2 = ndf.G(wi, wo, H);
            const float D = ndf.eval(H);

            // Ideal BRDF eval
            const float f = G2 * D / (4.f * wi.z);

            throughput += lobes[i].norm * f * factor * factor;
        }

        Assert(!throughput.isNaN());

        return throughput;
    }

    float pdf(const BSDFSamplingRecord& bRec, EMeasure measure) const {

        if (measure != ESolidAngle)
            return 0.f;

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        float       alphas[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, alphas);

        hg      lobes[TM_LAYER_COUNT_MAX];
        Vector3 lobes_wi[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, alphas, lobes, lobes_wi);

        ////////////////////////////////////////////////////////////////////////
        // PDF
        ////////////////////////////////////////////////////////////////////////

        float w_sum = 0.f;
        float wpdf_sum = 0.f;
        for (int i = 0; i < m_lobe_count; ++i) {

            if (lobes[i].norm.isZero())
                continue;

            const float norm = lobes[i].norm.average();

            w_sum += norm;

            if (lobes_wi[i].z * bRec.wo.z <= 0.f)
                continue;

            // Reflection eval directions (upper hemisphere)
            Vector3 wi = lobes_wi[i];
            Vector3 wo = bRec.wo;

            wi.z = std::abs(wi.z);
            wo.z = std::abs(wo.z);

            // Half-vector (upper hemisphere)
            const Vector3 H = normalize(wi + wo);

            // GGX equivalent normal distribution
            const float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const float G1  = ndf.smithG1(wi, H);
            const float D   = ndf.eval(H);

            wpdf_sum += norm * G1 * D / (4.0f * wi.z);
        }

        Assert(!(std::isnan(wpdf_sum) || std::isinf(wpdf_sum)));

        return w_sum > 0.f ? wpdf_sum / w_sum : 0.f;
    }

    Spectrum sample(BSDFSamplingRecord& bRec, const Point2& sample) const {

        float dummy;
        return this->sample(bRec, dummy, sample);
    }
    Spectrum sample(BSDFSamplingRecord& bRec, float& pdf, const Point2& sample) const {

        ////////////////////////////////////////////////////////////////////////
        // Lobes
        ////////////////////////////////////////////////////////////////////////

        // Get layers properties
        Spectrum    etas[TM_LAYER_COUNT_MAX];
        float       alphas[TM_LAYER_COUNT_MAX];

        evalMaps(bRec, etas, alphas);

        hg      lobes[TM_LAYER_COUNT_MAX];
        Vector3 lobes_wi[TM_LAYER_COUNT_MAX];

        outgoing_lobes(bRec.wi, etas, alphas, lobes, lobes_wi);

        //////////////////////////////////////////////////////////////////////
        // Random lobe selection
        //////////////////////////////////////////////////////////////////////

        // Lobes probabilities
        float w[TM_LAYER_COUNT_MAX];
        float w_sum = 0.f;
        for(int i = 0; i < m_lobe_count; ++i) {

            w[i] = lobes[i].norm.average();
            w_sum += w[i];
        }

        // Lobe selection
        float sel_w = bRec.sampler->next1D() * w_sum - w[0];
        int   sel_i = 0;
        for(sel_i = 0; sel_w > 0.f && sel_i < m_lobe_count; ++sel_i)
            sel_w -= w[sel_i + 1];

        //////////////////////////////////////////////////////////////////////
        // Sampling
        //////////////////////////////////////////////////////////////////////

        // Interaction type
        const bool reflection = bRec.wi.z * lobes_wi[sel_i].z > 0.f;

        // Selected lobe equivalent GGX normal distribution
        const float sel_a = hg_to_ggx(lobes[sel_i].g);
        const MicrofacetDistribution sel_ndf(MicrofacetDistribution::EGGX, sel_a, true);

        // Selected lobe reflection eval direction (upper hemisphere)
        Vector3 sel_wi = lobes_wi[sel_i];
        sel_wi.z = std::abs(sel_wi.z);

        // Selected lobe normal sampling (upper hemisphere)
        const Normal H = sel_ndf.sample(sel_wi, sample, pdf);

        // Sampled direction
        Vector3 sel_wo = reflect(sel_wi, H);
        sel_wo.z *= math::signum(bRec.wi.z) * (reflection ? 1.f : -1.f);

        // Sampling record
        bRec.wo = sel_wo;
        bRec.eta = etas[reflection ? 0 : m_layer_count].average();
        bRec.sampledComponent = reflection ? 0 : 1;
        bRec.sampledType = reflection ? EGlossyReflection : EGlossyTransmission;

        if(pdf <= 0.f)
            return Spectrum(0.0f);

        //////////////////////////////////////////////////////////////////////
        // PDF
        //////////////////////////////////////////////////////////////////////

        pdf = 0.f;
        for(int i = 0; i < m_lobe_count; ++i) {

            if (w[i] == 0.f)
                continue;

            if (lobes_wi[i].z * bRec.wo.z <= 0.f)
                continue;

            // Reflection eval directions (upper hemisphere)
            Vector3 wi = lobes_wi[i];
            Vector3 wo = bRec.wo;

            wi.z = std::abs(wi.z);
            wo.z = std::abs(wo.z);

            // Half-vector (upper hemisphere)
            const Vector3 H = normalize(wi + wo);

            // GGX equivalent normal distribution
            const float a = hg_to_ggx(lobes[i].g);
            const MicrofacetDistribution ndf(MicrofacetDistribution::EGGX, a, true);

            const float G1 = ndf.smithG1(wi, H);
            const float D = ndf.eval(H);

            pdf += w[i] * G1 * D / (4.0f * wi.z);
        }
        pdf /= w_sum;

        //////////////////////////////////////////////////////////////////////
        // Throughput
        //////////////////////////////////////////////////////////////////////

        Spectrum throughput = this->eval(bRec, ESolidAngle, etas, alphas);

        Assert(!throughput.isNaN());

        return pdf > 0.f ? throughput / pdf : Spectrum(0.0f);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(TM2_dielectric, false, BSDF)
MTS_EXPORT_PLUGIN(TM2_dielectric, "TM2_dielectric");
MTS_NAMESPACE_END

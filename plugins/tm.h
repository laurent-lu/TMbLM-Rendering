#ifndef _TM_
#define _TM_

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

#include "ior.h"
#include "microfacet.h"

#define TM_LAYER_COUNT_MAX 5
#define TM_EXP_ARG_MAX 45.f

#define CLAMP(v, vmin, vmax) ((std::max)((std::min)((v), (vmax)), (vmin)))

MTS_NAMESPACE_BEGIN

////////////////////////////////////////////////////////////////////////////////
// XML parsing and generic LUT class
//
// Code adapted from the supplemental code provided with "Efficient Rendering of
// Layered Materials using an Atomic Decomposition with Statistical Operators"
// Laurent Belcour SIGGRAPH 2018 paper
//
////////////////////////////////////////////////////////////////////////////////

void parseLayers(const Properties &props,
                 int& nb_layers,
                 std::vector<ref<const Texture>>& tex_etas,
                 std::vector<ref<const Texture>>& tex_kappas,
                 std::vector<ref<const Texture>>& tex_alphas,
                 bool opaque_material) {

    nb_layers = props.getInteger("nb_layers", 0);

    float extEta = lookupIOR(props, "extEta", "air");

    tex_etas.push_back(new ConstantSpectrumTexture(Spectrum(extEta)));
    tex_kappas.push_back(new ConstantSpectrumTexture(Spectrum(0.f)));
    tex_alphas.push_back(new ConstantSpectrumTexture(Spectrum(0.f)));

    for(int k=0; k<nb_layers; ++k) {

        std::string index = std::to_string(k);
        std::string name;

        name = std::string("eta_") + index;
        Spectrum eta = props.getSpectrum(name, Spectrum(extEta));

        name = std::string("kappa_") + index;
        Spectrum kappa = props.getSpectrum(name, Spectrum(0.f));

        // Check structure integrity
        if (!kappa.isZero()) {
            if (opaque_material) {
                if (k < (nb_layers - 1))
                    throw std::logic_error("Opaque material: Conductive layer allowed only as substrate");
            }
            else
                throw std::logic_error("Transparent material: Conductive layer not allowed");
        }
        else if (opaque_material && k == (nb_layers - 1))
            throw std::logic_error("Opaque material: Conductive substrate expected");

        name = std::string("alpha_") + index;
        float alpha = props.getFloat(name, 0.f);

        tex_etas.push_back(new ConstantSpectrumTexture(eta));
        tex_kappas.push_back(new ConstantSpectrumTexture(kappa));
        tex_alphas.push_back(new ConstantFloatTexture(alpha));

    }

    SLog(EInfo, "");
    SLog(EInfo, "Adding Exterior IOR");
    SLog(EInfo, " + n = %s", tex_etas[0]->getAverage().toString().c_str());
    SLog(EInfo, " + k = %s", tex_kappas[0]->getAverage().toString().c_str());
    SLog(EInfo, " + a = %f", tex_alphas[0]->getAverage().average());
    SLog(EInfo, "");

    for(int k=0; k<nb_layers; ++k) {

        SLog(EInfo, "Adding layer %d", k);
        SLog(EInfo, " + n = %s", tex_etas[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + k = %s", tex_kappas[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + a = %f", tex_alphas[k+1]->getAverage().average());
        SLog(EInfo, "");
    }
}
void parseLayers(const Properties &props,
                 int& nb_layers,
                 std::vector<ref<const Texture>>& tex_etas,
                 std::vector<ref<const Texture>>& tex_kappas,
                 std::vector<ref<const Texture>>& tex_alphas,
                 std::vector<ref<const Texture>>& tex_depths,
                 std::vector<ref<const Texture>>& tex_sigmas_s,
                 std::vector<ref<const Texture>>& tex_sigmas_k,
                 std::vector<ref<const Texture>>& tex_gs) {

    nb_layers = props.getInteger("nb_layers", 0);

    tex_etas.push_back(new ConstantSpectrumTexture(Spectrum(lookupIOR(props, "extEta", "air"))));
    tex_kappas.push_back(new ConstantSpectrumTexture(Spectrum(0.0)));
    tex_alphas.push_back(new ConstantFloatTexture(0.0));

    tex_depths.push_back(new ConstantFloatTexture(0.0));
    tex_sigmas_s.push_back(new ConstantSpectrumTexture(Spectrum(0.0)));
    tex_sigmas_k.push_back(new ConstantSpectrumTexture(Spectrum(0.0)));
    tex_gs.push_back(new ConstantFloatTexture(0.0));

    for(int k=0; k<nb_layers; ++k) {

        std::string index = std::to_string(k);
        std::string name;

        name = std::string("eta_") + index;
        Spectrum eta = props.getSpectrum(name, Spectrum(0.f));
        if (eta == Spectrum(0.f))
            eta = tex_etas.back()->getAverage();

        name = std::string("kappa_") + index;
        Spectrum kappa = props.getSpectrum(name, Spectrum(0.f));

        // Check structure integrity
        if (!kappa.isZero()) {
            if (k < (nb_layers - 1))
                throw std::logic_error("TM6: Conductive layer allowed only as substrate");
        }
        else if (k == (nb_layers - 1)) {
            throw std::logic_error("TM6: Conductive substrate expected");
        }

        name = std::string("alpha_") + index;
        float alpha = props.getFloat(name, 0.f);

        name = std::string("depth_") + index;
        float depth = props.getFloat(name, 0.f);

        name = std::string("sigmas_") + index;
        Spectrum sigma_s = props.getSpectrum(name, Spectrum(0.f));

        // Force non null absorption values
        name = std::string("sigmak_") + index;
        Spectrum sigma_k = props.getSpectrum(name, Spectrum(0.f));

        name = std::string("g_") + index;
        float g = props.getFloat(name, 1.f);

        tex_etas.push_back(new ConstantSpectrumTexture(eta));
        tex_kappas.push_back(new ConstantSpectrumTexture(kappa));
        tex_alphas.push_back(new ConstantFloatTexture(alpha));

        tex_depths.push_back(new ConstantFloatTexture(depth));
        tex_sigmas_s.push_back(new ConstantSpectrumTexture(sigma_s));
        tex_sigmas_k.push_back(new ConstantSpectrumTexture(sigma_k));
        tex_gs.push_back(new ConstantFloatTexture(g));
    }

    SLog(EInfo, "");
    SLog(EInfo, "Adding Exterior IOR");
    SLog(EInfo, " + n = %s", tex_etas[0]->getAverage().toString().c_str());
    SLog(EInfo, " + k = %s", tex_kappas[0]->getAverage().toString().c_str());
    SLog(EInfo, " + a = %f", tex_alphas[0]->getAverage().average());
    SLog(EInfo, "");

    for(int k=0; k<nb_layers; ++k) {

        SLog(EInfo, "Adding layer %d", k);

        SLog(EInfo, " + n = %s", tex_etas[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + k = %s", tex_kappas[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + a = %f", tex_alphas[k+1]->getAverage().average());
        SLog(EInfo, "");
        SLog(EInfo, " + d  = %f", tex_depths[k+1]->getAverage().average());
        SLog(EInfo, " + ss = %s", tex_sigmas_s[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + sk = %s", tex_sigmas_k[k+1]->getAverage().toString().c_str());
        SLog(EInfo, " + g  = %f", tex_gs[k+1]->getAverage().average());
        SLog(EInfo, "");
    }
}

struct lut_range {

	float min;
	float max;

	lut_range(float _min = std::numeric_limits<float>::quiet_NaN(),
			  float _max = std::numeric_limits<float>::quiet_NaN()) {

		min = _min;
		max = _max;
	}

	bool operator!=(const lut_range& other) const {

		return min != other.min || max != other.max;
	}
};

template<int LutDimension, typename LutDataType>
class lut {
private:

    // Data size
    int _size[LutDimension];
    // Data range
    lut_range _range[LutDimension];
    // Data
    std::vector<LutDataType> _data;

public:

    lut() { }
    lut(const std::string& filename) {

        std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);

        // Warn if not loaded
        if (in.bad() || in.fail())
            fprintf(stderr, "Unable to load %dD LUT data file: %s\n", LutDimension, filename.c_str());
        else
            fprintf(stdout, "Loading %dD LUT data file: %s\n", LutDimension, filename.c_str());

        this->load(in);
    }
    lut(const int size[], const lut_range range[]) {

        for (int i = 0; i < LutDimension; ++i) {

            _size[i] = size[i];
            _range[i] = range[i];
        }

        _data.assign(linear_size(), LutDataType());
    }

    const int* size() const { return _size; }
    constexpr int dimensions() const { return LutDimension; }
    const lut_range* range() const { return _range; }

    LutDataType* data() { return _data.data(); }
    const LutDataType* data() const { return _data.data(); }

    std::string size_string() const {

        std::ostringstream oss;

        for (int i = 0; i < LutDimension; ++i) {
            oss << _size[i];
            if (i < (LutDimension - 1))
                oss << "x";
        }

        return oss.str();
    }
    std::string range_string() const {

        std::ostringstream oss;
        for (int i = 0; i < LutDimension; ++i) {

            oss << "[" << _range[i].min << "," << _range[i].max << "]";
            if (i < (LutDimension - 1))
                oss << "x";
        }
        return oss.str();
    }

    int linear_size() const {

        int lsize = 1;
        for (int i = 0; i < LutDimension; ++i)
            lsize *= _size[i];
        return lsize;
    }
    int linear_index_from_grid_index(const int grid_index[]) const {

        int index = grid_index[LutDimension - 1];
        for (int i = LutDimension - 2; i >= 0; --i)
            index = index * _size[i] + grid_index[i];
        return index;
    }

    void load(std::ifstream& in) {

        // Read size
        for (int i = 0; i < LutDimension; ++i)
            in.read((char*)&_size[i], sizeof(int));

        fprintf(stdout, "Loading %dD LUT of dimensions %s\n", LutDimension, size_string().c_str());

        // Read data range (min / max)
        for (int i = 0; i < LutDimension; ++i) {

            in.read((char*)&_range[i].min, sizeof(float));
            in.read((char*)&_range[i].max, sizeof(float));
        }

        fprintf(stdout, "Loading %dD LUT data range: %s\n", LutDimension, range_string().c_str());

        // Read data
        _data.assign(linear_size(), LutDataType());
        in.read((char*)_data.data(), linear_size() * sizeof(LutDataType));
    }

    template<typename... T>
    inline LutDataType range_get(T... vrange_coords) const {

        float	x[LutDimension] = { vrange_coords... };
        int		i[LutDimension];

        for (int k = 0; k < LutDimension; ++k) {

            // Integer index (ensure the index stays in the limits)
            i[k] = (int)floor((_size[k] - 1) * (x[k] - _range[k].min) / (_range[k].max - _range[k].min));
            i[k] = CLAMP(i[k], 0, _size[k] - 1);
        }

        return _data[linear_index_from_grid_index(i)];
    }
    template<typename... T>
    inline LutDataType range_get_interpolate(T... vrange_coords) const {

        float	x[LutDimension] = { vrange_coords... };
        float	a[LutDimension];
        int		i[LutDimension];
        float	alphas[LutDimension];

        for (int k = 0; k < LutDimension; ++k) {

            // Floating point index
            a[k] = (_size[k] - 1) * (x[k] - _range[k].min) / (_range[k].max - _range[k].min);

            // Integer index
            i[k] = (int)floor(a[k]);

            // Ensure the indexes stays in the limits
            i[k] = CLAMP(i[k], 0, _size[k] - 1);

            // Clamp the interpolation weights
            alphas[k] = CLAMP(a[k] - (float)i[k], 0.f, 1.f);
        }

        int index = linear_index_from_grid_index(i);

        // Lookup result
        LutDataType v = LutDataType();

        // For every possible combinaison of index shift per dimension,
        // fetch the value in memory and do linear interpolation.
        // We fetch using shift of 0 and 1.
        //
        //     v(i+di, j+di, k+dk, l+dl),  where dk in [0,1]
        //
        const unsigned int D = (unsigned int)pow(2, LutDimension);
        for (unsigned int d = 0; d < D; ++d) {

            float alpha = 1.0; // Global alpha
            int   cid_s = 0;   // Id shift

            // Evaluate the weight of the sample d which correspond to
            // one for the shifted configuration:
            // The weight is the product of the weights per dimension.
            for (int k = (LutDimension - 1); k >= 0; --k) {

                bool  bitset = ((1 << k) & d);
                float calpha = (bitset) ? alphas[k] : 1.f - alphas[k];

                // Correct the shift to none if we go out of the grid
                if (i[k] + 1 >= _size[k])
                    bitset = false;

                alpha *= calpha;
                cid_s = cid_s * _size[k] + ((bitset) ? 1 : 0);
            }

            v += alpha * _data[index + cid_s];
        }

        return v;
    }
};

typedef lut<2, float> lut2;
typedef lut<3, float> lut3;
typedef lut<4, float> lut4;

////////////////////////////////////////////////////////////////////////////////
// Common routines
////////////////////////////////////////////////////////////////////////////////

static inline Vector3 reflectZ(const Vector3& v) {

    return Vector3(-v.x, -v.y, v.z);
}
static inline Vector3 refractZ(const Vector3& v, Float eta) {

    return refract(v, Vector3(0.f, 0.f,std::copysign(1.f, v.z)), eta);
}

static inline Spectrum min(const Spectrum& s, Float v) {

    Spectrum r;
    for (int i = 0; i < Spectrum::dim; ++i)
        r[i] = std::min(s[i], v);
    return r;
}
static inline Spectrum max(const Spectrum& s, Float v) {

    Spectrum r;
    for (int i = 0; i < Spectrum::dim; ++i)
        r[i] = std::max(s[i], v);
    return r;
}
static inline Spectrum exp(const Spectrum& s) {

    Spectrum r;
    for (int i = 0; i < Spectrum::dim; ++i)
        r[i] = std::exp(s[i]);
    return r;
}
static inline Spectrum sinh(const Spectrum& s) {

    Spectrum r;
    for (int i = 0; i < Spectrum::dim; ++i)
        r[i] = std::sinh(s[i]);
    return r;
}
static inline Spectrum safe_div(const Spectrum& n, const Spectrum& d) {

    Spectrum r;
    for (int i = 0; i < Spectrum::dim; ++i)
        r[i] = std::abs(d[i]) > 1e-7f ? n[i] / d[i] : 0.f;
    return r;
}

////////////////////////////////////////////////////////////////////////////////
// Henyey-Greenstein
////////////////////////////////////////////////////////////////////////////////

struct hg {

    // Asymmetry parameter
    Float g;
    // Mean
    Vector3 mean;
    // Norm
    Spectrum norm;

    hg(Float _g = 0.f, const Vector3& _mean = Vector3(0.f), const Spectrum& _norm = Spectrum(0.f))
        : g(_g), mean(_mean), norm(_norm) { }
};

static inline Float hg_from_ggx(Float a) {
    return CLAMP(-0.085f + (1.f + 0.085f)/(1.f + std::pow(std::min(CLAMP(a, 0.f, 1.f), 1.f)/0.5f, 1.3f)), 0.f, 1.f);
}
static inline Float hg_to_ggx(Float g) {
    return CLAMP(0.5f * std::pow((1.f + 0.085f)/(std::max(CLAMP(g, -1.f, 1.f), 0.2f) + 0.085f) - 1.f, 1.f / 1.3f), 1e-4f, 1.f);
}
static inline Float hg_refract(Float g, Float eta) {

  return std::min(std::sqrt(std::max(1.f - (1.f - g * g) * std::pow(eta, 0.75f), 0.f)), 1.f);
}
static inline Float hg_lh_norm(Float g) {

    const bool g_neg = g < 0.f;
    g = std::abs(g);
    const Float n = CLAMP(0.5039f - 0.8254f * g + 0.3226f * g * g, 0.f, 0.5f);
    return g_neg ? 1.f - n : n;
}

////////////////////////////////////////////////////////////////////////////////
// Albedos
////////////////////////////////////////////////////////////////////////////////

static inline void albedo(const Float cti, const Float alpha,
                          const Spectrum& etas_ij, const Spectrum& kappas_ij,
                          Spectrum& r_ij,
                          const lut4& FGD) {

    for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
        r_ij[i] = FGD.range_get_interpolate(cti, alpha, etas_ij[i], kappas_ij[i]);
}
static inline void albedos(const Float cti, const Float alpha, const Float eta_ij,
                           Spectrum& r_ij, Spectrum& t_ij, Spectrum& r_ji, Spectrum& t_ji,
                           const lut4& FGD) {

    // No IoR change
    if (std::abs(eta_ij - 1.f) < 1e-3f) {

        r_ij = r_ji = Spectrum(0.f);
        t_ij = t_ji = Spectrum(1.f);
        return;
    }

    r_ij = r_ji = Spectrum(FGD.range_get_interpolate(cti, alpha, eta_ij));
    t_ij = t_ji = Spectrum(1.f) - r_ij;
}

////////////////////////////////////////////////////////////////////////////////
// Common
////////////////////////////////////////////////////////////////////////////////

enum EComponentType : unsigned char {

    EUndefinedComponent,
    ENoComponent,
    EDielectricInterface,
    EConductorInterface,
    EHomogeneousMedium
};

////////////////////////////////////////////////////////////////////////////////
// 2-flux transfer calculus
////////////////////////////////////////////////////////////////////////////////

struct component_factors_2 {

    EComponentType type;

    hg r;     // Downward reflection
    hg t;     // Downward transmission
    hg rp;    // Upward reflection
    hg tp;    // Upward transmission

    component_factors_2(EComponentType type_ = EUndefinedComponent): type(type_) { }
};

template<typename Type>
struct transfer_matrix_2 {

    Type _m11, _m12;
    Type _m21, _m22;

    transfer_matrix_2():    _m11(Type(0.f)), _m12(Type(0.f)),
                            _m21(Type(0.f)), _m22(Type(0.f)) { }

    Type r() const { return _m21 / _m11; }
    Type t() const { return Type(1.f) / _m11; }
    
    Type r(const Type& r_cond) const {
        
        return (_m21 + _m22 * r_cond) / (_m11 + _m12 * r_cond);
    }

    inline transfer_matrix_2 operator*(const transfer_matrix_2& tm) const {

        transfer_matrix_2 res;

        res._m11 = _m11 * tm._m11 + _m12 * tm._m21;
        res._m12 = _m11 * tm._m12 + _m12 * tm._m22;
        
        res._m21 = _m21 * tm._m11 + _m22 * tm._m21;
        res._m22 = _m21 * tm._m12 + _m22 * tm._m22;

        return res;
    }
    inline transfer_matrix_2& operator*=(const transfer_matrix_2& tm) {

        *this = *this * tm;
        return *this;
    }

    static inline transfer_matrix_2 identity() {

        transfer_matrix_2 id;

        id._m11 = id._m22 = Type(1.f);

        return id;
    }
};

typedef transfer_matrix_2<Float>    transfer_matrix_21;
typedef transfer_matrix_2<Spectrum> transfer_matrix_2s;

static inline transfer_matrix_2s energy_matrix(const component_factors_2& ops) {

    transfer_matrix_2s etm;
    
    // Matrix factor
    const Spectrum t_inv = Spectrum(1.f) / ops.t.norm;

    // Matrix entries
    etm._m11 = t_inv;
    etm._m12 = -ops.rp.norm * t_inv;

    etm._m21 = ops.r.norm * t_inv;
    etm._m22 = (-ops.r.norm * ops.rp.norm + ops.t.norm * ops.tp.norm) * t_inv;

    return etm;
}
static inline transfer_matrix_21 asymmetry_matrix(const component_factors_2& ops) {

    transfer_matrix_21 gtm;

    // Weighted asymmetry parameters
    const Float r_g = ops.r.g * ops.r.norm.average();
    const Float t_g = ops.t.g * ops.t.norm.average();
    const Float rp_g = ops.rp.g * ops.rp.norm.average();
    const Float tp_g = ops.tp.g * ops.tp.norm.average();
    
    // Matrix factor
    const Float t_inv = 1.f / t_g;

    // Matrix entries
    gtm._m11 = t_inv;
    gtm._m12 = -rp_g * t_inv;

    gtm._m21 = r_g * t_inv;
    gtm._m22 = (-r_g * rp_g + t_g * tp_g) * t_inv;

    return gtm;
}

static inline void dielectric_transfer_factors(const Vector3& wl,
                                               const Float eta_ij,
                                               const Float alpha,
                                               const lut4& FGD,
                                               component_factors_2& ops) {

    Float eta_ji;
    Float s_t_ij, s_t_ji;

    // No IoR change
    if (std::abs(eta_ij - 1.f) < 1e-5f) {

        ops.type = ENoComponent;
        ops.t.mean = -wl;

        return;
    }

    // Component type
    ops.type = EDielectricInterface;

    // Inverse IoR
    eta_ji = 1.f / eta_ij;

    // Downward reflection
    ops.r.g = hg_from_ggx(alpha);
    ops.r.mean = reflectZ(wl);

    // Upward reflection
    ops.rp.g = ops.r.g;
    ops.rp.mean = refractZ(wl, eta_ij);

    // Fake rough transmissions scaling factors
    s_t_ij = std::abs((eta_ji * ops.r.mean.z + ops.rp.mean.z) / ops.rp.mean.z);
    s_t_ji = std::abs((eta_ij * ops.rp.mean.z + ops.r.mean.z) / ops.r.mean.z);

    // Downward transmission (fake reflection)
    ops.t.g = hg_from_ggx(0.5f * s_t_ij * alpha);
    ops.t.mean = ops.rp.mean;

    // Upward transmission (fake reflection)
    ops.tp.g = hg_from_ggx(0.5f * s_t_ji * alpha);
    ops.tp.mean = ops.r.mean;

    // Albedos
    albedos(std::abs(wl.z), alpha, eta_ij, ops.r.norm, ops.t.norm, ops.rp.norm, ops.tp.norm, FGD);
}
static inline void conductor_transfer_factors(const Vector3& wl,
                                              const Spectrum& etas,
                                              const Spectrum& kappas,
                                              const Float alpha,
                                              const lut4& FGD,
                                              component_factors_2& ops) {

    // Component type
    ops.type = EConductorInterface;

    // Downward reflection
    ops.r.g = hg_from_ggx(alpha);
    ops.r.mean = reflectZ(wl);

    // Albedo
    albedo(std::abs(wl.z), alpha, etas, kappas, ops.r.norm, FGD);
}

////////////////////////////////////////////////////////////////////////////////
// 6-flux transfer calculus
////////////////////////////////////////////////////////////////////////////////

struct component_factors_6 {

    EComponentType type;

    union {

        // Interface transfer factors
        struct {

            hg r;       // Primary downward reflection
            hg t;       // Primary downward transmission
            hg rp;      // Primary upward reflection
            hg tp;      // Primary upward transmission

            hg r_ff;    // Secondary downward reflection
            hg t_ff;    // Secondary downward transmission
            hg rp_ff;   // Secondary upward reflection
            hg tp_ff;   // Secondary upward transmission
        } i;
        // Homogeneous participating medium transfer factors
        struct {

            hg t;     // Primary flux transmission
            hg r_fb;  // Secondary flux backward reflection
            hg t_ff;  // Secondary flux forward transmission
        } m;
    };

    component_factors_6(EComponentType type_ = EUndefinedComponent): type(type_) { }
};

enum ETransferFactorRequest : unsigned char {

    ENone                           = 0x00,
    EPrimaryReflectance             = 0x01,
    ESecondaryForwardReflectance    = 0x02,
    ESecondaryBackwardReflectance   = 0x04
};

template<typename Type>
struct transfer_matrix_6 {

    EComponentType type;

    Type m11, m12;
    Type m21, m22;
    Type m31, m32, m33, m34, m35, m36;
    Type m41, m42, m43, m44, m45, m46;
    Type m51, m52;
    Type m61, m62;

    transfer_matrix_6()
     : type(EUndefinedComponent),
       m11(Type(0.f)), m12(Type(0.f)),
       m21(Type(0.f)), m22(Type(0.f)),
       m31(Type(0.f)), m32(Type(0.f)), m33(Type(0.f)), m34(Type(0.f)), m35(Type(0.f)), m36(Type(0.f)),
       m41(Type(0.f)), m42(Type(0.f)), m43(Type(0.f)), m44(Type(0.f)), m45(Type(0.f)), m46(Type(0.f)),
       m51(Type(0.f)), m52(Type(0.f)),
       m61(Type(0.f)), m62(Type(0.f)) { }

    static inline transfer_matrix_6 identity() {

        transfer_matrix_6 id;

        id.type = ENoComponent;
        id.m11 = id.m22 = id.m33 = id.m44 = Type(1.f);

        return id;
    }

    inline transfer_matrix_6 operator*(const transfer_matrix_6& tm) const {

        transfer_matrix_6 res;

        // Identity
        if (type == ENoComponent)
            return tm;

        if (tm.type == ENoComponent) {

            return *this;
        }
        else if (tm.type == EDielectricInterface) {

            res.m11 = m11*tm.m11 + m12*tm.m21;
            res.m12 = m11*tm.m12 + m12*tm.m22;

            res.m21 = m21*tm.m11 + m22*tm.m21;
            res.m22 = m21*tm.m12 + m22*tm.m22;

            res.m31 = m31*tm.m11 + m32*tm.m21;
            res.m32 = m31*tm.m12 + m32*tm.m22;
            res.m33 = m33*tm.m33 + m34*tm.m43;
            res.m34 = m33*tm.m34 + m34*tm.m44;
            res.m35 = m35*tm.m33 + m36*tm.m43;
            res.m36 = m35*tm.m34 + m36*tm.m44;

            res.m41 = m41*tm.m11 + m42*tm.m21;
            res.m42 = m41*tm.m12 + m42*tm.m22;
            res.m43 = m43*tm.m33 + m44*tm.m43;
            res.m44 = m43*tm.m34 + m44*tm.m44;
            res.m45 = m45*tm.m33 + m46*tm.m43;
            res.m46 = m45*tm.m34 + m46*tm.m44;

            res.m51 = m51*tm.m11 + m52*tm.m21;
            res.m52 = m51*tm.m12 + m52*tm.m22;

            res.m61 = m61*tm.m11 + m62*tm.m21;
            res.m62 = m61*tm.m12 + m62*tm.m22;
        }
        else if (tm.type == EHomogeneousMedium) {

            res.m11 = m11*tm.m11;
            res.m12 = m12*tm.m22;

            res.m21 = m21*tm.m11;
            res.m22 = m22*tm.m22;

            res.m31 = m31*tm.m11 + m33*tm.m31 - m36*tm.m36;
            res.m32 = m32*tm.m22 + m34*tm.m42 + m35*tm.m36;
            res.m33 = m33*tm.m33 - m36*tm.m36;
            res.m34 = m34*tm.m44 + m35*tm.m36;
            res.m35 = -m34*tm.m36 + m35*tm.m33;
            res.m36 = m33*tm.m36 + m36*tm.m44;

            res.m41 = m41*tm.m11 + m43*tm.m31 - m46*tm.m36;
            res.m42 = m42*tm.m22 + m44*tm.m42 + m45*tm.m36;
            res.m43 = m43*tm.m33 - m46*tm.m36;
            res.m44 = m44*tm.m44 + m45*tm.m36;
            res.m45 = -m44*tm.m36 + m45*tm.m33;
            res.m46 = m43*tm.m36 + m46*tm.m44;

            res.m51 = -m34*tm.m36 + m35*tm.m31 + m51*tm.m11;
            res.m52 = m33*tm.m36 + m36*tm.m42 + m52*tm.m22;

            res.m61 = -m44*tm.m36 + m45*tm.m31 + m61*tm.m11;
            res.m62 = m43*tm.m36 + m46*tm.m42 + m62*tm.m22;
        }
        else {
            throw std::logic_error("Operation not supported");
        }

        return res;
    }
    inline transfer_matrix_6& operator*=(const transfer_matrix_6& tm) {

        *this = *this * tm;
        return *this;
    }

    inline void factors(Type& r, Type& r_f, Type& r_b, unsigned char req) const {

        const Type x0 = Type(1.f)/m11;

        if (req & ESecondaryBackwardReflectance) {

            const Type x1 = -m31*m33 + m35*m51;
            const Type x2 = m31*m35 - m33*m51;
            const Type x3 = x0/(m33*m33 - m35*m35);
            const Type x4 = m43*x3;
            const Type x5 = m45*x3;

            r_b = m61*x0 + x1*x5 + x4*x2;
            if (req & ESecondaryForwardReflectance)
                r_f = m41*x0 + x1*x4 + x2*x5;
        }

        if (req & EPrimaryReflectance)
            r = m21*x0;
    }
    inline void factors(Type& r, Type& r_f, Type& r_b, const Type& r_cond, unsigned char req) const {

        const Type x0 = Type(1.f)/(m11 + m12*r_cond);

        if (req & ESecondaryBackwardReflectance) {

            const Type x1 = m45 + m46*r_cond;
            const Type x2 = m31 + m32*r_cond;
            const Type x3 = m35 + m36*r_cond;
            const Type x4 = m33 + m34*r_cond;
            const Type x5 = m51 + m52*r_cond;
            const Type x6 = x0/(-x3*x3 + x4*x4);
            const Type x7 = x6*(x2*x3 - x4*x5);
            const Type x8 = m43 + m44*r_cond;
            const Type x9 = x6*(x2*x4 - x3*x5);

            r_b = x0*(m61 + m62*r_cond) - x1*x9 + x7*x8;
            if (req & ESecondaryForwardReflectance)
                r_f = x0*(m41 + m42*r_cond) + x1*x7 - x8*x9;
        }

        if (req & EPrimaryReflectance)
            r = x0*(m21 + m22*r_cond);
    }
};

typedef transfer_matrix_6<Spectrum> transfer_matrix_6s;

static inline transfer_matrix_6s energy_matrix(const component_factors_6& ops) {

    transfer_matrix_6s etm;

    etm.type = ops.type;

    if (ops.type == EDielectricInterface) {

        etm.m11 = Spectrum(1.f) / ops.i.t.norm;
        etm.m12 = -ops.i.rp.norm * etm.m11;
        etm.m21 = ops.i.r.norm * etm.m11;
        etm.m22 = -ops.i.r.norm * ops.i.rp.norm * etm.m11 + ops.i.tp.norm;

        etm.m33 = Spectrum(1.f) / ops.i.t_ff.norm;
        etm.m34 = -ops.i.rp_ff.norm * etm.m33;
        etm.m43 = ops.i.r_ff.norm * etm.m33;
        etm.m44 = -ops.i.r_ff.norm * ops.i.rp_ff.norm * etm.m33 + ops.i.tp_ff.norm;
    }
    else if (ops.type == EHomogeneousMedium) {

        etm.m11 = Spectrum(1.f) / ops.m.t.norm;
        etm.m22 = ops.m.t.norm;
        etm.m33 = Spectrum(1.f) / ops.m.t_ff.norm;
        etm.m44 = -(ops.m.r_fb.norm * ops.m.r_fb.norm) * etm.m33 + ops.m.t_ff.norm;
        etm.m52 = etm.m36 = -ops.m.r_fb.norm * etm.m33;
        etm.m31 = etm.m33 - etm.m11;
        etm.m42 = etm.m44 - ops.m.t.norm;
        etm.m61 = etm.m45 = -etm.m36;
    }
    else {
        throw std::logic_error("Operation not supported");
    }

    return etm;
}
static inline transfer_matrix_6s asymmetry_matrix(const component_factors_6& ops) {

    transfer_matrix_6s gtm;

    gtm.type = ops.type;

    if (ops.type == EDielectricInterface) {

        const Spectrum r_g     = ops.i.r.norm  * ops.i.r.g;
        const Spectrum t_g     = ops.i.t.norm  * ops.i.t.g;
        const Spectrum rp_g    = ops.i.rp.norm * ops.i.rp.g;
        const Spectrum tp_g    = ops.i.tp.norm * ops.i.tp.g;

        const Spectrum r_ff_g   = ops.i.r_ff.norm  * ops.i.r_ff.g;
        const Spectrum t_ff_g   = ops.i.t_ff.norm  * ops.i.t_ff.g;
        const Spectrum rp_ff_g  = ops.i.rp_ff.norm * ops.i.rp_ff.g;
        const Spectrum tp_ff_g  = ops.i.tp_ff.norm * ops.i.tp_ff.g;

        gtm.m11 = Spectrum(1.f) / t_g;
        gtm.m12 = -rp_g * gtm.m11;
        gtm.m21 = r_g * gtm.m11;
        gtm.m22 = -r_g * rp_g * gtm.m11 + tp_g;

        gtm.m33 = Spectrum(1.f) / t_ff_g;
        gtm.m34 = -rp_ff_g * gtm.m33;
        gtm.m43 = r_ff_g * gtm.m33;
        gtm.m44 = -r_ff_g * rp_ff_g * gtm.m33 + tp_ff_g;
    }
    else if (ops.type == EHomogeneousMedium) {

        const Spectrum m_t_g       = ops.m.t.norm    * ops.m.t.g;
        const Spectrum m_r_fb_g    = ops.m.r_fb.norm * ops.m.r_fb.g;
        const Spectrum m_t_ff_g    = ops.m.t_ff.norm * ops.m.t_ff.g;

        gtm.m11 = Spectrum(1.f) / m_t_g;
        gtm.m22 = m_t_g;
        gtm.m33 = Spectrum(1.f) / m_t_ff_g;
        gtm.m44 = -(m_r_fb_g * m_r_fb_g) * gtm.m33 + m_t_ff_g;
        gtm.m52 = gtm.m36 = -m_r_fb_g * gtm.m33;
        gtm.m31 = gtm.m33 - gtm.m11;
        gtm.m42 = gtm.m44 - m_t_g;
        gtm.m61 = gtm.m45 = -gtm.m36;
    }
    else {
        throw std::logic_error("Operation not supported");
    }

    return gtm;
}

static inline void dielectric_transfer_factors(const Vector3& wl,
                                               const Float eta_ij,
                                               const Float alpha,
                                               const lut4& FGD,
                                               component_factors_6& ops) {

	Float eta_ji;
    Float s_t_ij, s_t_ji;

    // No IoR change
    if (std::abs(eta_ij - 1.f) < 1e-5f) {

        ops.type = ENoComponent;
        ops.i.t.mean = -wl;

        return;
    }

    // Component type
    ops.type = EDielectricInterface;

    // Inverse IoR
    eta_ji = 1.f / eta_ij;

    // Downward reflection
    ops.i.r.g = hg_from_ggx(alpha);
    ops.i.r.mean = reflectZ(wl);

    // Upward reflection
    ops.i.rp.g = ops.i.r.g;
    ops.i.rp.mean = refractZ(wl, eta_ij);

    // Fake rough transmissions scaling factors
    s_t_ij = std::abs((eta_ji * ops.i.r.mean.z + ops.i.rp.mean.z) / ops.i.rp.mean.z);
    s_t_ji = std::abs((eta_ij * ops.i.rp.mean.z + ops.i.r.mean.z) / ops.i.r.mean.z);

    // Downward transmission (fake reflection)
    ops.i.t.g = hg_from_ggx(0.5f * s_t_ij * alpha);
    ops.i.t.mean = ops.i.rp.mean;

    // Upward transmission (fake reflection)
    ops.i.tp.g = hg_from_ggx(0.5f * s_t_ji * alpha);
    ops.i.tp.mean = ops.i.r.mean;

    // Albedos
    albedos(std::abs(wl.z), alpha, eta_ij, ops.i.r.norm, ops.i.t.norm, ops.i.rp.norm, ops.i.tp.norm, FGD);

    ops.i.r_ff     = ops.i.r;
    ops.i.t_ff     = ops.i.t;
    ops.i.rp_ff    = ops.i.rp;
    ops.i.tp_ff    = ops.i.tp;
}
static inline void conductor_transfer_factors(const Vector3& wl,
                                              const Spectrum& etas,
                                              const Spectrum& kappas,
                                              const Float alpha,
                                              const lut4& FGD,
                                              component_factors_6& ops) {

    // Component type
    ops.type = EConductorInterface;

    // Downward reflection
    ops.i.r.g = hg_from_ggx(alpha);
    ops.i.r.mean = reflectZ(wl);

    // Albedo
    albedo(std::abs(wl.z), alpha, etas, kappas, ops.i.r.norm, FGD);

    ops.i.r_ff = ops.i.r;
}
static inline void medium_transfer_factors(const Vector3& wl,
                                           const Float depth,
                                           const Spectrum& sigma_s,
                                           const Spectrum& sigma_k,
                                           const Float g,
                                           component_factors_6& ops) {

    // Component type
    ops.type = EHomogeneousMedium;

    // Empty medium
    if (depth == 0.f) {

        ops.type = ENoComponent;
        ops.i.t.mean = -wl;

        return;
    }

    // Incident depth
    Float tau = depth / wl.z;

    // Interaction coefficients
    const Spectrum sigma_sb     = sigma_s * hg_lh_norm(g);
    const Spectrum sigma_sf     = sigma_s - sigma_sb;
    const Spectrum sigma_ext    = sigma_s + sigma_k;

    // Energies
    {
        const Spectrum alpha    = sigma_ext - sigma_sf;
        const Spectrum beta     = sigma_sb;
        const Spectrum gamma    = (alpha * alpha - beta * beta).sqrt();

        // Clamp depth to avoid numerical issues with transcendent functions
        tau = std::min(tau, TM_EXP_ARG_MAX / gamma.max());
        tau = std::min(tau, TM_EXP_ARG_MAX / sigma_ext.max());

        const Spectrum S = sinh(gamma * tau);
        const Spectrum C = (Spectrum(1.f) + S * S).sqrt();

        const Spectrum x0 = exp(-sigma_ext * tau);
        const Spectrum x1 = Spectrum(1.f) / (C * gamma + S * alpha);

        ops.m.t.norm     = x0;
        ops.m.r_fb.norm  = S * beta * x1;
        ops.m.t_ff.norm  = gamma * x1;
    }
    // Asymetry parameters
    {
        const Spectrum alpha   = sigma_ext - sigma_sf * g;
        const Spectrum beta    = sigma_sb * -g;
        const Spectrum gamma   = (alpha * alpha - beta * beta).sqrt();

        const Spectrum S = sinh(gamma * tau);
        const Spectrum C = (Spectrum(1.f) + S * S).sqrt();

        const Spectrum x0 = Spectrum(1.f) / (C * gamma + S * alpha);

        ops.m.t.g    = 1.f;
        ops.m.r_fb.g = safe_div(S * beta * x0, ops.m.r_fb.norm).average();
        ops.m.t_ff.g = safe_div(gamma * x0, ops.m.t_ff.norm).average();
    }
}

MTS_NAMESPACE_END

#endif

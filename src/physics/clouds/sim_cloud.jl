export SimCloud
"""
Cloud fraction is a linear function of relative humidity.
A freeze-dry adjustment is applied to the cloud fraction to compensate for 
excessive cloudiness over the poles, especially the winter pole.
A low cloud fraction adjustment is applied to the cloud fraction to compensate for
underestimation of marine low cloud.

References:
- SimCloud version 1.0: a simple diagnostic cloud scheme for idealized climate models
  Liu, Q., Collins, M., Maher, P., Thomson, S. I., Vallis, G. K.
  Geoscientific Model Development, 14(5), 2801-2826, 2021
  DOI: https://dio.org/10.5194/gmd-14-2801-2021
  URL: https://gmd.copernicus.org/articles/14/2801/2021/

$(TYPEDFIELDS)"""
@kwdef struct SimCloud{NF<:AbstractFloat} <: AbstractLargeScaleCloud
    "a at surface"
    a_s::NF = 36.0 # [dimensionless]

    "a in the free troposphere"
    a_t::NF = 13.0 # [dimensionless]

    "n used in calculation of the gradient of the cf-rh curve"
    n::Int = 12 # [dimensionless]

    "the power used to describe how quickly specific humidity drops with height in the freeze-dry adjustment"
    n_freeze_dry::NF = 2.5 # [dimensionless]

    "the proscribed surface specific humidity used to control the freeze-dry adjustment"
    q0::NF = 0.006 # [kg/kg]

    "the proscribed sea level pressure used to control the freeze-dry adjustment"
    pres_msl::NF = 100000.0 # [Pa]

    "the threshold lapse rate used to control the low cloud fraction adjustment"
    dthetadp_threshold::NF = -0.00125 # [K/Pa]

    "Cloud cover scheme"
    cloud_cover_scheme::CloudCover = RandomOverlap()
end

SimCloud(SG::SpectralGrid; kwargs...) = SimCloud{SG.NF}(; kwargs...)
function initialize!(scheme::SimCloud, model::PrimitiveEquation)
    # initialise cloud cover scheme
    initialize!(scheme.cloud_cover_scheme, model)
end

function large_scale_cloud!( 
    column::ColumnVariables,
    scheme::SimCloud,
    model::PrimitiveWet,
)
    large_scale_cloud!(column, scheme)
    cloud_cover!(column, scheme.cloud_cover_scheme)
end

function large_scale_cloud!(
    column::ColumnVariables,
    scheme::SimCloud,
)
    (; pres, humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; bulk_cloud_fraction) = column        # diagnostics to write into
    (; sat_humid, temp_virt) = column       # intermediate variables, calculated in thermodynamics!
    
    (; a_s, a_t, n, n_freeze_dry, q0, pres_msl, dthetadp_threshold) = scheme

    @inbounds for k in eachindex(column)
        # compute the gradient of the cloud fraction - relative humidity curve
        a = a_t + (a_s - a_t) * exp(1 - (pres_surf / pres[k])^n)
        # compute relative humidity
        RH = humid[k]/sat_humid[k]
        # compute cloud fraction
        bulk_cloud_fraction[k] = min(1, max(0, a * (RH - 1) + 1))
        # freeze-dry adjustment
        qv = q0 * (pres[k] / pres_msl)^n_freeze_dry
        bulk_cloud_fraction[k] *= max(0.15, min(1.0, (humid[k] / qv) ))
        # low cloud fraction adjustment (CURRENTLY SKIPPED)
        # if below 750hPa, get lapse rate and vertical velocity
        # if pres[k] < 75000
        #     # find most negative lapse rate dthetadp
        #     # TODO: check this before use: gradients in the right place, correct temperature used?
        #     # are we at the top of the model?
        #     if k == length(column)
        #         # set lapse rate to gradient between here and the layer below
        #         dthetadp = min(dthetadp, (temp_virt[k] - temp_virt[k-1]) / (pres[k] - pres[k-1]))
        #     else
        #         dthetadp = min(dthetadp, (temp_virt[k+1] - temp_virt[k]) / (pres[k+1] - pres[k]))
        #     end
        #     # get vertical velocity
        #     # TODO: vertical velocity ought to be in the column variables. From here it is assumed to be > 0 Pa/s.
        # end

        # If there is now clouds at a level higher (i.e. smaller k) than
        # the cloud-top previously diagnosed due to convection, then increase the cloud-top
        if bulk_cloud_fraction[k] > 0
            column.cloud_top = min(column.cloud_top, k)
        end
    end
end

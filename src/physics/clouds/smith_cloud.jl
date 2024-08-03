export SmithCloud
"""
Cloud cover assuming a symmetric triangular distribution of subgrid saturation departure.
We specify a critical relative humidity at the surface and top of the boundary layer, and a boundary layer height.
The critical relative humidity is linearly interpolated between these two points.

The original Smith cloud scheme did not directly diagnose convective clouds, nor does this implementation.
A subsequent improvement was the Empirically Adjusted Cloud Fraction (EACF) version of this scheme,
which boosted cf for profiles near saturation to better match observations. This is not implemented here.

References:
- Smith, R. N. B., 1990: A scheme for predicting layer clouds and their water content in a general circulation model. 
Quart. J. Roy. Meteor. Soc., 116, 435-460, https://doi.org/10.1002/qj.49711649210. 

- Weverberg, K.V., Morcrette, C.J., Boutle, I., Furtado, K. & Field, P.R. (2021):
A bimodal diagnostic cloud fraction parameterization. Part I: motivating analysis and scheme description. 
Monthly Weather Review, 149, 841-857, https://doi.org/10.1175/MWR-D-20-0224.1. 
 
$(TYPEDFIELDS)"""
@kwdef struct SmithCloud{NF<:AbstractFloat} <: AbstractLargeScaleCloud
    "RHcrit at surface"
    RHcrit_surf::NF = 0.96

    "RHcrit at top of boundary layer"
    RHcrit_top::NF = 0.8

    "Approximate boundary layer height [Pa]"
    BLH::NF = 90000.0

    "Cloud cover scheme"
    cloud_cover_scheme::CloudCover = RandomOverlap()
end

SmithCloud(SG::SpectralGrid; kwargs...) = SmithCloud{SG.NF}(; kwargs...)
function initialize!(scheme::SmithCloud, model::PrimitiveEquation)
    # initialise cloud cover scheme
    initialize!(scheme.cloud_cover_scheme, model)
end

function large_scale_cloud!( 
    column::ColumnVariables,
    scheme::SmithCloud,
    model::PrimitiveWet,
)
    large_scale_cloud!(column, scheme)
    cloud_cover!(column, scheme.cloud_cover_scheme)
end

function large_scale_cloud!(
    column::ColumnVariables,
    scheme::SmithCloud,
)
    (; pres, humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; bulk_cloud_fraction) = column        # diagnostics to write into
    (; sat_humid) = column                  # intermediate variable, calculated in thermodynamics!
    
    (; RHcrit_surf, RHcrit_top, BLH) = scheme

    @inbounds for k in eachindex(column)
        # linearly drop RHcrit from surface to top of boundary layer
        RHcrit = RHcrit_surf - (RHcrit_surf - RHcrit_top) * (pres[k] - pres[end]) / BLH
        # compute relative humidity wrt total water content (all vapour rn because speedy doesn't track other phases)
        RH = humid[k]/sat_humid[k]
        # cloud cover is 0 if RH <= RHcrit, 1 if 2 - RHcrit <= RH, 
        # Given by 1/2 * (1+(RH-1)/(1-RHcrit))^2 if RHcrit < RH <= 1,
        # and 1 - 1/2 * (1+(RH-RHcrit)/(1-RHcrit))^2 if 1 < RH < 2 - RHcrit
        if RH <= RHcrit
            bulk_cloud_fraction[k] = 0
        elseif RH >= 2 - RHcrit
            bulk_cloud_fraction[k] = 1
        elseif RH <= 1
            bulk_cloud_fraction[k] = 1/2 * (1+(RH-1)/(1-RHcrit))^2
        else
            bulk_cloud_fraction[k] = 1 - 1/2 * (1+(RH-RHcrit)/(1-RHcrit))^2
        end
        # If there is now clouds at a level higher (i.e. smaller k) than
        # the cloud-top previously diagnosed due to convection, then increase the cloud-top
        if bulk_cloud_fraction[k] > 0
            column.cloud_top = min(column.cloud_top, k)
        end
    end
end

export HomogeneousCloud
"""
Cloud cover present if grid box mean relative humidity is above 1.
$(TYPEDFIELDS)"""
@kwdef struct HomogeneousCloud{NF<:AbstractFloat} <: AbstractLargeScaleCloud
    "Relative humidity threshold [1 = 100%] to trigger cloud cover"
    relative_humidity_threshold::NF = 1

    "cloud cover scheme (default: MaximumOverlap as bulk cloud fraction in this scheme is either 0 or 1)"
    cloud_cover_scheme::CloudCover = MaximumOverlap()
end
HomogeneousCloud(SG::SpectralGrid; kwargs...) = HomogeneousCloud{SG.NF}(; kwargs...)
function initialize!(scheme::HomogeneousCloud, model::PrimitiveEquation)
    # initialise cloud cover scheme
    initialize!(scheme.cloud_cover_scheme, model)
end

function large_scale_cloud!( 
    column::ColumnVariables,
    scheme::HomogeneousCloud,
    model::PrimitiveWet,
)
    large_scale_cloud!(column, scheme)
    cloud_cover!(column, scheme.cloud_cover_scheme)
    
end

function large_scale_cloud!(
    column::ColumnVariables,
    scheme::HomogeneousCloud,
)
    (; humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; bulk_cloud_fraction) = column        # diagnostics to write into
    (; sat_humid) = column                  # intermediate variable, calculated in thermodynamics!
    
    @inbounds for k in eachindex(column)
        # cloud cover is 0 if RH is below threshold, else 1.
        if humid[k]/sat_humid[k] < scheme.relative_humidity_threshold
            bulk_cloud_fraction[k] = 0
        else
            bulk_cloud_fraction[k] = 1
        end
        # If there is now clouds at a level higher (i.e. smaller k) than
        # the cloud-top previously diagnosed due to convection, then increase the cloud-top
        if bulk_cloud_fraction[k] > 0
            column.cloud_top = min(column.cloud_top, k)
        end
    end
end

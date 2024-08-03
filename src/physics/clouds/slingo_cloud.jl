export SlingoCloud
"""
Cloud cover as a quadratic function of relative humidity. 
Clouds are divided into three classes based on altitude.

References:
- Slingo, J. M. (1987). The development and verification of a cloud prediction scheme for the ECMWF model. 
  Quarterly Journal of the Royal Meteorological Society, 113(477), 899-927.
$(TYPEDFIELDS)"""
@kwdef struct SlingoCloud{NF<:AbstractFloat} <: AbstractLargeScaleCloud
    "RH threshold for low level cloud (1=100%)"
    Mlow::NF = 0.8
    
    "RH threshold for mid level cloud (1=100%)"
    Mmid::NF = 0.65
    
    "RH threshold for high level cloud (1=100%)"
    Mhigh::NF = 0.8
    
    "Pressure threshold for low level cloud (Pa)"
    plow::NF = 80000.0

    "Pressure threshold for high level cloud (Pa)"
    phigh::NF = 40000.0

    "Cloud cover scheme"
    cloud_cover_scheme::CloudCover = RandomOverlap()
end    
SlingoCloud(SG::SpectralGrid; kwargs...) = SlingoCloud{SG.NF}(; kwargs...)
function initialize!(scheme::SlingoCloud, model::PrimitiveEquation)
    print("hey i am the slingo cloud scheme edited with cloud cover")
    # initialise cloud cover scheme
    initialize!(scheme.cloud_cover_scheme, model)
end


# function barrier for SlingoCloud to unpack model
function large_scale_cloud!( 
    column::ColumnVariables,
    scheme::SlingoCloud,
    model::PrimitiveWet,
)
    large_scale_cloud!(column, scheme)
    cloud_cover!(column, scheme.cloud_cover_scheme)
end

"""
$(TYPEDSIGNATURES)
Diagnostic large-scale cloud for a `column` by assuming the saturation
departure follows the subgrid distribution ((RH-M)/(1-M))^2 where M is a constant"""
function large_scale_cloud!(
    column::ColumnVariables,
    scheme::SlingoCloud,
)
    (; pres, humid) = column          # prognostic vars (from previous time step for numerical stability)
    (; bulk_cloud_fraction) = column        # diagnostics to write into
    (; sat_humid) = column                  # intermediate variable, calculated in thermodynamics!
    
    (; Mlow, Mmid, Mhigh, plow, phigh) = scheme

    @inbounds for k in eachindex(column)
        # check which of our three cloud types the pressure here falls into
        p = pres[k]
        if p > plow
            M = Mlow
        elseif p < phigh
            M = Mhigh
        else
            M = Mmid
        end
        # cloud cover is 0 if RH is below M, 1 if RH>=1, and (RH-M)/(1-M)² in between
        RH = humid[k]/sat_humid[k]
        if RH < M
            bulk_cloud_fraction[k] = 0
        elseif RH >= 1
            bulk_cloud_fraction[k] = 1
        else
            bulk_cloud_fraction[k] = ((RH-M)/(1-M))^2
        end
        # TODO: the original Slingo scheme diagnoses convective cloud cover based on convective precipitation
        #       and convective cloud base and top. This is not implemented here as we currently don't have
        #       cloud base diagnosed by any convection schemes. Thus convective clouds are ignored.
        #       The convective cloud cover would affect both high and mid level clouds.
        #        (check what we can do with only cloud top)
        
        # TODO: the original Slingo scheme splits low clouds into two types based on vertical velocity,
        #       with a small transition between them.
        #       This is not implemented here as we currently don't have vertical velocity in the ColumnVariables.
        #       Thus this scheme simply assumes omega >= 0.
        # If there is now clouds at a level higher (i.e. smaller k) than
        # the cloud-top previously diagnosed due to convection, then increase the cloud-top
        if bulk_cloud_fraction[k] > 0
            column.cloud_top = min(column.cloud_top, k)
        end
    end
end


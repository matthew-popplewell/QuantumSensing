#=
Synthetic star fields calibrated to Gaia DR3 star counts.

Star counts per 1-mag G bin were measured once from Gaia DR3 in 1-deg cones
on a grid of galactic (l, b) (see data/gaia_dr3_star_counts.csv for the
query). Synthetic stars are placed at random positions with numbers drawn
from those counts (interpolated to the pointing). Stars fainter than the
drawn limit are summed into an integrated-starlight surface brightness, so
the resolved stars and the diffuse starlight stay consistent.
=#

const STAR_COUNT_FILE = joinpath(@__DIR__, "..", "data", "gaia_dr3_star_counts.csv")
const STAR_COUNT_CONE_DEG2 = 2pi * (1 - cosd(1.0)) * (180 / pi)^2  # 1-deg cone
const STAR_G_BINS = collect(0:27)      # bin k covers G in [k, k + 1)
const STAR_MAX_SLOPE = 0.38            # cap on faint-end d log10(N) / dG

"""
    SyntheticStar

A synthetic star: `ra`, `dec` (Radians, J2000) and V magnitude `v_mag`.
"""
struct SyntheticStar
    ra::Radians
    dec::Radians
    v_mag::Float64
end

# Lazily loaded count grid: Dict((l, b) => (density per deg^2, mean BP-RP)).
const _STAR_COUNTS = Ref{Any}(nothing)

function _load_star_counts()
    _STAR_COUNTS[] === nothing || return _STAR_COUNTS[]
    fields = Dict{Tuple{Float64, Float64}, Tuple{Vector{Float64}, Vector{Float64}}}()
    for line in eachline(STAR_COUNT_FILE)
        (startswith(line, "#") || startswith(line, "l_deg")) && continue
        l, b, g, n, c = parse.(Float64, split(line, ","))
        dens, col = get!(fields, (l, b)) do
            (zeros(length(STAR_G_BINS)), fill(NaN, length(STAR_G_BINS)))
        end
        k = Int(g) + 1
        dens[k] = n / STAR_COUNT_CONE_DEG2
        col[k] = c
    end
    for (dens, col) in values(fields)
        _extrapolate_faint!(dens, col)
    end
    _STAR_COUNTS[] = fields
    return fields
end

"""
    _extrapolate_faint!(density, colour)

Replace Gaia-incomplete faint bins with an exponential extrapolation (as Krag
2003 did): the counts are trusted up to the last bin before they stop
rising, then continued with the mean log-slope of the three preceding bins
(capped at `STAR_MAX_SLOPE`). The colour of the last trusted bin is reused.
"""
function _extrapolate_faint!(density, colour)
    first_faint = findfirst(==(15), STAR_G_BINS)
    last_ok = first_faint
    while last_ok < length(density) - 1 && density[last_ok + 1] > density[last_ok]
        last_ok += 1
    end
    last_ok = min(last_ok, findfirst(==(20), STAR_G_BINS))
    slope = (log10(density[last_ok]) - log10(density[last_ok - 3])) / 3
    slope = clamp(slope, 0.0, STAR_MAX_SLOPE)
    for k in (last_ok + 1):length(density)
        density[k] = density[last_ok] * 10.0^(slope * (k - last_ok))
        colour[k] = colour[last_ok]
    end
    # Bins with no stars have no colour; use a typical value (it is unused
    # because their density is zero).
    replace!(colour, NaN => 1.0)
    return density
end

"""
    gaia_g_minus_v(bp_rp)

Gaia DR3 G - V (Johnson-Cousins) from BP - RP, valid for -0.5 < BP-RP < 5
(sigma 0.030 mag). Gaia DR3 documentation, Sec. 5.5.1, Table 5.9:
https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html
"""
function gaia_g_minus_v(bp_rp)
    c = clamp(bp_rp, -0.5, 5.0)
    return -0.02704 + 0.01424 * c - 0.2156 * c^2 + 0.01426 * c^3
end

"""
    galactic_coordinates(direction)

Galactic longitude and latitude (degrees) of a J2000 unit vector.
"""
function galactic_coordinates(direction)
    g = pxform("J2000", "GALACTIC", 0.0) * collect(direction)
    return mod(atand(g[2], g[1]), 360.0), asind(clamp(g[3] / norm(g), -1.0, 1.0))
end

"""
    star_counts(direction)

Star surface density per 1-mag G bin toward a J2000 unit `direction`.

# Returns
NamedTuple with `g` (lower edge of each bin), `density` (stars per deg^2 in
the bin) and `v_offset` (V - G for the bin's mean colour).
"""
function star_counts(direction)
    fields = _load_star_counts()
    l, b = galactic_coordinates(direction)
    lats = sort(unique(last.(collect(keys(fields)))))
    lons = sort(unique(first.(k for k in keys(fields) if abs(k[2]) < 90)))

    # Log-density and colour of one latitude row, interpolated in longitude.
    function row(lat)
        abs(lat) == 90 && return fields[(0.0, lat)]
        dl = 360.0 / length(lons)
        i = clamp(Int(fld(l, dl)) + 1, 1, length(lons))
        w = clamp((l - lons[i]) / dl, 0.0, 1.0)
        l2 = lons[mod1(i + 1, length(lons))]
        d1, c1 = fields[(lons[i], lat)]
        d2, c2 = fields[(l2, lat)]
        return (_log_mix(d1, d2, w), (1 - w) .* c1 .+ w .* c2)
    end

    j, w = _bracket(lats, b)
    d1, c1 = row(lats[j])
    d2, c2 = row(lats[j + 1])
    density = _log_mix(d1, d2, w)
    colour = (1 - w) .* c1 .+ w .* c2
    return (g = STAR_G_BINS, density = density, v_offset = -gaia_g_minus_v.(colour))
end

# Interpolate densities in log space (zeros stay zero only if both are zero).
function _log_mix(d1, d2, w)
    return [
        (a > 0 && c > 0) ? 10.0^((1 - w) * log10(a) + w * log10(c)) : (1 - w) * a + w * c
        for (a, c) in zip(d1, d2)
    ]
end

"""
    integrated_starlight_s10(counts, v_limit)

Surface brightness (S10, above the atmosphere) of all stars in bins whose
mean V is at or fainter than `v_limit`.
"""
function integrated_starlight_s10(counts, v_limit::Real)
    s10 = 0.0
    for (g, n, dv) in zip(counts.g, counts.density, counts.v_offset)
        v = g + 0.5 + dv
        v >= v_limit && (s10 += n * 10.0^(-0.4 * (v - 10)))
    end
    return s10
end

"""
    synthetic_star_field(center, radius, v_limit; rng=Random.default_rng())

Random stars in a cone of angular `radius` around the J2000 unit vector
`center`, for all count bins whose mean V is brighter than `v_limit`. Counts
per bin are Poisson-distributed; magnitudes are uniform within each bin.
"""
function synthetic_star_field(
    center::AbstractVector,
    radius,
    v_limit::Real;
    rng::AbstractRNG = Random.default_rng()
)
    rho = ustrip(u"rad", to_unit(u"rad", radius))
    c = normalize(SVector{3, Float64}(center))
    counts = star_counts(c)
    area_deg2 = 2pi * (1 - cos(rho)) * (180 / pi)^2

    # Orthonormal basis around the field centre.
    helper = abs(c[3]) < 0.9 ? SVector(0.0, 0.0, 1.0) : SVector(1.0, 0.0, 0.0)
    e1 = normalize(cross(helper, c))
    e2 = cross(c, e1)

    stars = SyntheticStar[]
    for (g, n, dv) in zip(counts.g, counts.density, counts.v_offset)
        (n > 0 && g + 0.5 + dv < v_limit) || continue
        for _ in 1:rand(rng, Poisson(n * area_deg2))
            cos_t = 1 - rand(rng) * (1 - cos(rho))
            sin_t = sqrt(1 - cos_t^2)
            phi = 2pi * rand(rng)
            u = cos_t * c + sin_t * (cos(phi) * e1 + sin(phi) * e2)
            push!(stars, SyntheticStar(
                atan(u[2], u[1]) * u"rad",
                asin(clamp(u[3], -1.0, 1.0)) * u"rad",
                g + rand(rng) + dv
            ))
        end
    end
    return stars, counts
end

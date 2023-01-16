# **************************************************************************** #
#                                                                              #
#                          Temperature functions                               #
#                                                                              #
# **************************************************************************** #

function T(Tsurf::Float64, Tmeso::Float64, Texo::Float64; z_meso_top=110e5, lapserate=-8e-5, weird_Tn_param=8, globvars...)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    
    Uses the typical temperature structure for neutrals, but interpolates temperatures
    for ions and electrons above 108 km according to the profiles in Fox & Sung 2001.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:alt])
    
    # Subroutines -------------------------------------------------------------------------------------

    function upperatmo_i_or_e(z, particle_type; Ti_array=Ti_interped, Te_array=Te_interped, select_alts=new_a)
        #=
        Finds the index for altitude z within new_a
        =#
        i = something(findfirst(isequal(z), select_alts), 0)
        returnme = Dict("electron"=>Te_array[i], "ion"=>Ti_array[i])
        return returnme[particle_type]
    end

    function NEUTRALS()
        function upper_atmo_neutrals(z_arr)
            @. return Texo - (Texo - Tmeso)*exp(-((z_arr - z_meso_top)^2)/(weird_Tn_param*1e10*Texo))
        end
        Tn = zeros(size(GV.alt))

        Tn[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Tn[i_meso] .= Tmeso
        Tn[i_upper] .= upper_atmo_neutrals(GV.alt[i_upper])

        return Tn 
    end 

    function ELECTRONS(;spc="electron") 
        Te = zeros(size(GV.alt))

        Te[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Te[i_meso] .= Tmeso

        # Upper atmo
        Telec = readdlm("../Resources/PlotsFromPapers/foxandsung2001_Te.dat",' ', Float64, comments=true, comment_char='#')
        T_e = Telec[:, 1]
        alt_e = Telec[:, 2] .* 1e5
        interp_elec = LinearInterpolation(alt_e, T_e)
        Te_interped = [interp_elec(a) for a in new_a];
        Te[i_upper] .= Te_interped# upperatmo_i_or_e(alt[i_upper], spc)

        return Te
    end

    function IONS(;spc="ion") 
        Ti = zeros(size(GV.alt))

        Ti[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Ti[i_meso] .= Tmeso

        # Upper atmo
        Tion = readdlm("../Resources/PlotsFromPapers/foxandsung2001_Tion.dat",' ', Float64, comments=true, comment_char='#')
        T_i = Tion[:, 1]
        alt_i = Tion[:, 2] .* 1e5
        interp_ion = LinearInterpolation(alt_i, T_i)
        Ti_interped = [interp_ion(a) for a in new_a];
        Ti[i_upper] .= Ti_interped# upperatmo_i_or_e(alt[i_upper], spc)

        return Ti
    end

    # Define mesosphere scope.
    z_meso_bottom = alt[searchsortednearest(alt, (Tmeso-Tsurf)/(lapserate))]

    # Various indices to define lower atmo, mesosphere, and atmo -----------------
    i_lower = findall(z->z < z_meso_bottom, GV.alt)
    i_meso = findall(z->z_meso_bottom <= z <= z_meso_top, GV.alt)
    i_upper = findall(z->z > z_meso_top, GV.alt)
    # i_meso_top = findfirst(z->z==z_meso_top, GV.alt)

    # For interpolating upper atmo temperatures from Fox & Sung 2001
    new_a = collect(112e5:2e5:250e5)


    # In the lower atmosphere, neutrals, ions, and electrons all have the same temperatures. 
    # if z < z_meso_bottom
    #     return Tsurf + lapserate*z
    # elseif z_meso_bottom <= z <= z_meso_top 
    #     return Tmeso
    # # Near the top of the isothermal mesosphere, profiles diverge.        
    # elseif z > z_meso_top
    #     if sptype=="neutral"
    #         return T_upper_atmo_neutrals(z)
    #     else
    #         return upperatmo_i_or_e(z, sptype)
    #     end
    # end

    return Dict("neutrals"=>NEUTRALS(), "ions"=>IONS(), "electrons"=>ELECTRONS())
end
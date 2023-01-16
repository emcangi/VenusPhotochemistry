################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Global constants, simulation parameters.
# 
# Eryn Cangi
# Created December 2019
# Last edited: October 2022
# Currently tested for Julia: 1.7.1
################################################################################

using DataFrames

# **************************************************************************** #
#                                                                              #
#                         Main simulation parameters                           #
#                                                                              #
# **************************************************************************** #

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                                                                        !! #
# !!                      !!!!! SUPER IMPORTANT !!!!!                       !! #
# !!     !!! Modify the following items each time you run the model !!!     !! #
# !!                                                                        !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

reinitialize_atmo = true

# Basic simulation parameters
const optional_logging_note = "First Venus Attempt" # Simulation goal
const simset = "VenusPaper"
const results_version = "v0"  # Helps keep track of attempts 
const initial_atm_file = "INITIAL_GUESS.h5" 

const seasonal_cycle = false # False for Venus
const timestep_type =  "dynamic-log" #"log-linear" # basically never use this one: "static-log" 

# SET EXPERIMENT
const exptype =  "temperature" #"water" # "insolation"#
# temperature
# const use_for_Texo = 225.
#water
const water_case = "1e-6" # for Venus, this is just the mixing ratio.
const water_mixing_ratio = 1e-6 #
# const update_water_profile = false # this is for modifying the profile during cycling, currently doesn't work well
# solar cycle
const solarcyc = "mean" # "max" # "min"#

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                                                                        !! #
# !!                      !!!!!    END CHECK    !!!!!                       !! #
# !!                                                                        !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# Set up temperature and folder tag name according to the type of simulation set

const solarfile = "VENUSsolarphotonflux_solarmean.dat" # CHANGE
const controltemps = [735., 170., 260.]
const tag = "venus_temp_$(results_version)"

# Temperature and water

const meantemps = [735., 170., 260.] # Used for saturation vapor pressure. DON'T CHANGE!
const reinitialize_water_profile = true # false # setting to true will set the water profile = water_mixing_ratio everywhere

# Tolerance and timespans 
const season_length_in_sec = seasonal_cycle==true ? sol_in_sec : 1e16
const maxlogdt = seasonal_cycle==true ? 5 : 16
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-3, maxlogdt]) 
const rel_tol = 1e-6
const abs_tol = 1e-12 

# More basics that don't frequently change
const ions_included = true # false # Turn off when starting from neutrals
const SZA = 60 # degrees 
const fixed_species = [] # here you may enter any species that you want to be completely fixed (no updates to densities from chemistry or transport)

# Tags, shortcodes, and filenames
const hrshortcode, rshortcode = generate_code(ions_included, controltemps[1], controltemps[2], controltemps[3], water_case, solarcyc)
const sim_folder_name = "$(hrshortcode)_$(rshortcode)_$(tag)"
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"REACTION_NETWORK.xlsx" # "REACTION_NETWORK_MIN_IONOSPHERE.xlsx" #

# Ionospheric chemistry and non-thermal escape
const nontherm = ions_included==true ? true : false  
const converge_which = "both"
const e_profile_type = "none" # "quasineutral" #"O2+" # "constant"# Set to "none" when attempting to converge from neutrals first, quasineutral otherwise
const remove_unimportant = true # Whether to use a slightly smaller list of species and reactions (removing minor species that Roger had in his model)
const adding_new_species = false # set to true if introducing a new species.


# **************************************************************************** #
#                                                                              #
#                             Algorithm settings                               #
#                                                                              #
# **************************************************************************** #

# Algorithm and timestepping scheme
const problem_type = "Gear" #"SS" #"ODE" #  
# for static timesteps:
const n_steps = 800 # for static case
# for dynamic timesteps:
const dt_incr_factor = 1.5
const dt_decr_factor = 10
# other solver details:
const ediff = false # true # whether to include differentiation terms in Jacobian with respect to electron density
const mdiff = false # true # whether to include differentiation terms in Jacobian with respect to density of the generic thirdbody "M"
const error_checking_scheme = "new" #"old" 

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be
if problem_type == "Gear"
    const assume_photochem_eq = false
else
    const assume_photochem_eq = converge_which == "both" ? true : false
end

# Check that float type is properly set
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in CUSTOMIZATIONS.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end


# **************************************************************************** #
#                                                                              #
#      Misc. things used when adding new species or changing altitude grid     #
#                                                                              #
# **************************************************************************** #
const do_chem = true  # false # set to false to try diffusion only
const do_trans = true 
const make_new_alt_grid = false
const use_nonzero_initial_profiles = true

# **************************************************************************** #
#                                                                              #
#               Temperature profile construction from controls                 #
#                                                                              #
# **************************************************************************** #

const T_surf = controltemps[1]
const T_meso = controltemps[2]
const T_exo = controltemps[3]

T_array_dict = T(T_surf, T_meso, T_exo; alt);
const Tn_arr = T_array_dict["neutrals"]
const Ti_arr = T_array_dict["ions"]
const Te_arr = T_array_dict["electrons"]

const Tplasma_arr = Ti_arr .+ Te_arr;
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)
const Tn_meanSVP = T(meantemps...; alt)["neutrals"]; # Needed for boundary conditions.

# **************************************************************************** #
#                                                                              #
#                             Boundary conditions                              #
#                                                                              #
# **************************************************************************** #
const H2Osat = map(x->Psat(x), Tn_meanSVP) # Using this function keeps SVP fixed 
# const HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))

const KoverH_lowerbdy = Keddy([zmin], [ntot_at_lowerbdy])[1]/scaleH_lowerboundary(zmin, Tn_arr[1]; molmass)

const speciesbclist=Dict(:CO2=>Dict("n"=>[0.965*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0.]),
                         :CO=>Dict("n"=>[4.5e-6*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0.]),
                         :O=> Dict("v"=>[KoverH_lowerbdy, NaN], "f"=>[NaN, 1.2e6]),
                         :O2=>Dict("v"=>[KoverH_lowerbdy, NaN], "f"=>[NaN, 0]),
                         :N=>Dict("v"=>[KoverH_lowerbdy, NaN], "f"=>[NaN, 0]),
                         :NO=>Dict("n"=>[5.5e-9*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0]),

                         :H=> Dict("v"=>[NaN, effusion_velocity(Tn_arr[end], 1.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                         :D=> Dict("v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                         :H2=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),  # velocities are in cm/s
                         :HD=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 3.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),

                         :Hpl=>Dict("f"=>[NaN, 1.6e7]),
                         :H2pl=>Dict("f"=>[NaN, 2e5]),
                         :Opl=>Dict("f"=>[NaN, 5.2e6]),
                       );

# **************************************************************************** #
#                                                                              #
#                    Species name lists and J rate lists                       #
#                                                                              #
# **************************************************************************** #

unimportant = [:CNpl,:HCNpl,:HCNHpl,:HN2Opl,:NH2pl,:NH3pl,:N2Opl,:NO2pl,:CH,:CN,:HCN,:HNO,:NH,:NH2,:N2O,:NO2,:HD2pl]

# Neutrals --------------------------------------------------------------------
const orig_neutrals = [:CO2, :Ar, :CO, :CO2, :H, :H2, :H2O, :H2O2, 
                       :HO2, :HOCO, :N2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD,

                       # Turn these off for minimal ionosphere:
                       :C, :CH, :CN, :DCO, :HCN, :HCO, :HNO, :N, :NH, :NH2, :N2O, :NO, :NO2, :Nup2D,
                       ]; 
const conv_neutrals = remove_unimportant==true ? setdiff(orig_neutrals, unimportant) : orig_neutrals
const new_neutrals = [];
const neutral_species = [conv_neutrals..., new_neutrals...];

# Ions -------------------------------------------------------------------------
const orig_ions = [:CO2pl, :HCO2pl, :Opl, :O2pl, # Nair minimal ionosphere 
                   :Arpl, :ArHpl, :ArDpl, 
                   :Cpl, :CHpl, :CNpl, :COpl, 
                   :Hpl, :Dpl, :H2pl, :HDpl, :H3pl, :H2Dpl, :HD2pl, 
                   :H2Opl,  :HDOpl, :H3Opl, :H2DOpl, 
                   :HO2pl, :HCOpl, :DCOpl, :HOCpl, :DOCpl, :DCO2pl, 
                   :HCNpl, :HCNHpl, :HNOpl, :HN2Opl,  
                   :Npl,  :NHpl, :NH2pl, :NH3pl, :N2pl, :N2Hpl, :N2Dpl, :N2Opl, :NOpl, :NO2pl,
                   :OHpl, :ODpl];
const new_ions = []; 
# const ion_species = [conv_ions..., new_ions...];
const ion_species = remove_unimportant==true ? setdiff([orig_ions..., new_ions...], unimportant) : [orig_ions..., new_ions...]

# Full species list -------------------------------------------------------------
const all_species = [neutral_species..., ion_species...];

# Sort name lists created here -------------------------------------------------
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)

# Photolysis and Photoionization rate symbol lists ----------------------------

const conv_Jrates, newJrates = format_Jrates(reaction_network_spreadsheet, all_species, "Jratelist"; ions_on=ions_included, hot_atoms=nontherm)
const Jratelist = [conv_Jrates..., newJrates...];

# These dictionaries specify the species absorbing a photon for each J rate, and the products of the reaction.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# **************************************************************************** #
#                                                                              #
#                      Miscellaneous logical groupings                         #
#                                                                              #
# **************************************************************************** #
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = get_deuterated(all_species)
const D_ions = get_deuterated(ion_species) #[s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];

# Sort name lists created here -------------------------------------------------
sort!(D_bearing_species)
sort!(D_ions)
sort!(N_neutrals)

# **************************************************************************** #
#                                                                              #
#                    Define short- and long-lived species                      #
#                                                                              #
# **************************************************************************** #

# Short lived species, whose chemical lifetime is << diffusion timescale -------
const short_lived_species = [];
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2, :C, :CH])
    append!(short_lived_species, ion_species)
end

# Long lived species ------------------------------------------------------------
const long_lived_species = setdiff(all_species, short_lived_species)

# Sort name lists created here -------------------------------------------------
sort!(short_lived_species)
sort!(long_lived_species)

# **************************************************************************** #
#                                                                              #
#               Species participating in chemistry and transport               #
#                                                                              #
# **************************************************************************** #

# Non-participants -------------------------------------------------------------
const no_chem_species = []; 
const no_transport_species = [];

# This will append any species that you have picked to be completely fixed. 
for fs in fixed_species
    push!(no_chem_species, fs)
    push!(no_transport_species, fs)
end

if converge_which == "neutrals"
    append!(no_chem_species, union(conv_ions, N_neutrals)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_neutrals))
    append!(no_transport_species, setdiff(conv_neutrals, N_neutrals))
elseif converge_which == "both"
    append!(no_transport_species, short_lived_species)
end

# Participants ------------------------------------------------------------------
const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species ---------------------------------------------------
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# Sort name lists created here -------------------------------------------------
sort!(active_species)
sort!(inactive_species)
sort!(active_longlived)
sort!(active_shortlived)
sort!(chem_species)
sort!(transport_species)
sort!(no_chem_species)
sort!(no_transport_species)

# **************************************************************************** #
#                                                                              #
#              Misc. things that depend on things defined above                #
#                                                                              #
# **************************************************************************** #

# Annoyingly, these have to be here because they depend on other things defined above.
# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict(vcat([s=>"--" for s in setdiff(D_bearing_species, [:HD2pl])], [:HD2pl=>":", :Nup2D=>"-."]) )

# Species-specific scale heights - has to be done here instead of in the param file
const Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)]; molmass) for sp in all_species])

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position with the active species vector 
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

# Altitude at which water transitions from fixed to freely solved for
# First we have to calculate a few intermediaries.
H2Osatfrac = H2Osat ./ map(z->n_tot(get_ncurrent(initial_atm_file), z; all_species, n_alt_index), alt)  # get SVP as fraction of total atmo
const upper_lower_bdy = alt[something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)] # in cm
const upper_lower_bdy_i = n_alt_index[upper_lower_bdy]  # the uppermost layer at which water will be fixed, in cm

# **************************************************************************** #
#                                                                              #
#               CREATE A PARAMETER DATAFRAME FOR LOGGING EASE                  #
#                                                                              #
# **************************************************************************** # 

PARAMETERS_GEN = DataFrame(Field=[], Value=[])

push!(PARAMETERS_GEN, ("HRSHORTCODE", hrshortcode));
push!(PARAMETERS_GEN, ("RSHORTCODE", rshortcode));
push!(PARAMETERS_GEN, ("INITIAL_ATM", initial_atm_file));
push!(PARAMETERS_GEN, ("RXN_SOURCE", reaction_network_spreadsheet));
push!(PARAMETERS_GEN, ("IONS", ions_included ));
push!(PARAMETERS_GEN, ("CONVERGE", converge_which));
push!(PARAMETERS_GEN, ("NONTHERMAL_ESC", nontherm));
push!(PARAMETERS_GEN, ("SOLARCYCLE", solarcyc));
push!(PARAMETERS_GEN, ("SOLARFILE", solarfile));
push!(PARAMETERS_GEN, ("ELECTRON_PROF", e_profile_type));
push!(PARAMETERS_GEN, ("EDIFF", ediff));
push!(PARAMETERS_GEN, ("MDIFF", mdiff));
push!(PARAMETERS_GEN, ("DH", DH));

PARAMETERS_CONDITIONS = DataFrame(Field=[], Value=[], Unit=[]);

push!(PARAMETERS_CONDITIONS, ("SZA", SZA, "deg"));
push!(PARAMETERS_CONDITIONS, ("TSURF", T_surf, "K"));
push!(PARAMETERS_CONDITIONS, ("TMESO", T_meso, "K"));
push!(PARAMETERS_CONDITIONS, ("TEXO", T_exo, "K"));
push!(PARAMETERS_CONDITIONS, ("MEAN_TEMPS", join(meantemps, " "), "K"));
push!(PARAMETERS_CONDITIONS, ("WATER_MR", water_mixing_ratio, "mixing ratio"));
push!(PARAMETERS_CONDITIONS, ("WATER_CASE", water_case, "whether running with 10x, 1/10th, or standard water in middle/upper atmo"));
push!(PARAMETERS_CONDITIONS, ("WATER_BDY", upper_lower_bdy/1e5, "km"))

# This is so ugly because the XLSX package won't write columns of different lengths, so I have to pad all the shorter lists
# with blanks up to the length of the longest list and also transform all the symbols into strings. 
L = max(length(all_species), length(neutral_species), length(ion_species), length(no_chem_species), length(no_transport_species))
PARAMETERS_SPLISTS = DataFrame(AllSpecies=[[string(a) for a in all_species]..., ["" for i in 1:L-length(all_species)]...], 
                               Neutrals=[[string(n) for n in neutral_species]..., ["" for i in 1:L-length(neutral_species)]...], 
                               Ions=[[string(i) for i in ion_species]..., ["" for i in 1:L-length(ion_species)]...],
                               NoChem=[[string(nc) for nc in no_chem_species]..., ["" for i in 1:L-length(no_chem_species)]...],
                               NoTransport=[[string(nt) for nt in no_transport_species]..., ["" for i in 1:L-length(no_transport_species)]...],
                               Jratelist=[[string(j) for j in Jratelist]..., ["" for i in 1:L-length(Jratelist)]...]);

PARAMETERS_SOLVER = DataFrame(Field=[], Value=[]);
PARAMETERS_XSECTS = DataFrame(Species=[], Description=[], Filename=[]);
PARAMETERS_BCS = DataFrame(Species=[], Type=[], Lower=[], Upper=[]);


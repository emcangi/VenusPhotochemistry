# **************************************************************************** #
#                                                                              #
#                          Basic utility functions                             #
#                                                                              #
# **************************************************************************** #

# Checked for planet-neutrality 20-Dec-22

#                          Standard, miscellaneous functions                   #
#==============================================================================#

function deletefirst(A, v)
    #=
    returns: list A with its first element equal to v removed. Used to make the derivative for 
    the chemical jacboian 
    =#
    index = something(findfirst(isequal(v), A), 0)  # this horrible syntax introduced by Julia devs
    keep = setdiff([1:length(A);],index)
    return A[keep]
end

function df_lookup(df, indcol, indcolentry, col)
    #=
    df: Dataframe to search, produced by final_escape.
    indcol: The name of the index column, such as "EscapeType"
    indcolentry: The value you want in the index column, i.e. which row to use, such as 'Thermal' or 'Total'
    col: Column name, probably a species like "H"
    =#
    return df[in([indcolentry]).(df.:($indcol)), col]
end

function find_nonfinites(collection; collec_name="collection")
    #=
    Returns indices of any nonfinite values (inf or nan) in collection.
    =#
    nonfinites = findall(x->x==0, map(el->isfinite(el), collection))
    if length(nonfinites) != 0
        throw("ALERT: Found nonfinite values in $(collec_name) at indices $(nonfinites)")
    end
end

function fluxsymbol(x)
    #= 
    Converts string x to a symbol. f in the string means flux. 
    =#
    return Symbol(string("f",string(x)))
end

function generate_code(ii, TS, TM, TE, water, scyc; 
                       print_optional=false, pt="Gear", tstype="D", elecval="q", rt=1e-2, at=1e-12)
    #=
    Generates a relatively human-readable short code for identifying main simulation features
    and also returns a unique randomly generated shortcode for embeddingin plots and h5 files to
    identify the source simulation of a file.

    Input: These all represent global variables in the parameters file.
        pt: problem_type
        ii: ions_included
        tstype: timestep_type
        elecval: e_profile_type
        nt: nontherm
        TS, TM, TE: T_surf, T_meso, T_exo
        water: water state (high, low)
        scyc: solarcyc
        rt: rel_tol
        at: abs_tol
    =#
    iontag = ii == true ? "I" : "N"
    temptag = "Ts$(Int64(TS))Tm$(Int64(TM))Te$(Int64(TE))"
    watertag = "W$(water)"
    solartag = "S$(scyc)"
    
    tags = [iontag, temptag, watertag, solartag]

    if print_optional
        # optional tags
        timesteptag = Dict("static"=>"S", "dynamic"=>"D")[tstype]
        etag = Dict("constant"=>"cc", "O2+"=>"O2+", "quasineutral"=>"q", "none"=>"0")[elecval]
        rttag = "RT$(rt)"
        attag = "AT$(at)"
        esctag = nt == true ? "NT" : "T" 
        tags = [tags..., "_", timesteptag, pt, etag, rttag, attag, nt, esctag]
    end 

    hrcode = join(tags, "-")
    
    return hrcode, randstring()
end

function get_paramfile(working_dir; use_default=false)
    #=
    Helper function to let the user select a parameter file form working_dir. 
    =#
    available_paramfiles = filter(x->(occursin(".jl", x) & occursin("PARAMETERS", x)), readdir(working_dir))

    if length(available_paramfiles) == 1
        paramfile = available_paramfiles[1]
        println("Using the only file available, $(paramfile)")
    else 
        if use_default == true
            paramfile = "PARAMETERS.jl"
        else
            println("Available parameter files: ")
            [println("[$(i)] $(available_paramfiles[i])") for i in 1:length(available_paramfiles)]
            paramfile_selection = input("Select a parameter file or enter filename with extension: ")
            paramfile = available_paramfiles[parse(Int64, paramfile_selection)]
        end
        println("Using parameter file $(paramfile).")
    end

    return paramfile
end

function getpos(array, test::Function, n=Any[])
    #= 
    Get the position of keys that match "test" in an array

    array: Array; any size and dimension
    test: used to test the elements in the array
    n: Array; for storing the indices (I think)

    returns: 1D array of elements in array that match test function supplied
    =#
    if !isa(array, Array)
        test(array) ? Any[n] : Any[]
    else
        vcat([ getpos(array[i], test, Any[n...,i]) for i=1:size(array)[1] ]...)
    end
end

function getpos(array, value)
    #= 
    overloading getpos for the most common use case, finding indicies
    corresponding to elements in array that match value. 
    =#
    getpos(array, x->x==value)
end

function input(prompt::String="")::String
    #=
    Prints a prompt to the terminal and accepts user input, returning it as a string.
    =#
    print(prompt)
    return chomp(readline())
end

function logrange(x1, x2, n::Int64)
    #=
    Equivalent to python's logspace, returns a generator with entries that
    are powers of 10.
    x1: start, some power of 10
    x2: end, some power of 10
    n: number of entries, spaced between the powers of x1 and x2

    e.g.: collect(logrange(1e-4, 1e2, 4)) returns 0.0001, 0.01, 1.0, 100.0

    =#
    return (10.0^y for y in range(log10(x1), log10(x2), length=n))
end 

nans_present(a) = any(x->isnan(x), a)

function next_in_loop(i::Int64, n::Int64)
    #=
    returns i+1, restricted to values within [1, n]. 
    next_in_loop(n, n) returns 1.
    =#
    return i % n + 1
end

function searchsortednearest(a,x)
    #=
    This does something. Mike wrote it. I think it's used in the cross sections.
    =#
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end

function subtract_difflength(a::Array, b::Array)
    #=
    A very specialized function that accepts two vectors, a and b, sorted
    by value (largest first), of differing lengths. It will subtract b from a
    elementwise up to the last index where they are equal, and then add any 
    extra values in a, and subtract any extra values in b.

    Used exclusively in ratefn_local. 

    a: production rates 
    b: loss rates

    both a and b should be positive for the signs to work!
    =#

    shared_size = min(size(a), size(b))[1]

    extra_a = 0
    extra_b = 0
    if shared_size < length(a)
        extra_a += sum(a[shared_size+1:end])
    end

    if shared_size < length(b)
        extra_b += sum(b[shared_size+1:end])
    end

    return sum(a[1:shared_size] .- b[1:shared_size]) + extra_a - extra_b
end

#                        String manipulation functions                          #
#===============================================================================#

function charge_type(sp::Symbol)
    #=
    Returns a string representing the type of sp, i.e. ion, neutral, or electron
    =#
    if occursin("pl", String(sp))
        return "ion"
    elseif sp==:E
        return "electron"
    else
        return "neutral"
    end
end

function decompose_chemistry_string(s; returntype="array")
    #=
    Takes a formatted chemistry string, such as that produced by format_chemistry_string, 
    and returns an array of vectors of the form [reactants, products].

    e.g. "OH + OH --> H2O + O" returns [[:OH, :OH], [:H2O, :O]]
    =#
    reactants, products = split(s, " --> ")
    if returntype=="strings"
        return reactants, products
    end

    reactants = [Symbol(s) for s in split(reactants, " + ")]
    products = [Symbol(s) for s in split(products, " + ")]

    if returntype=="array"
        return [reactants, products]
    else
        throw("Bad returntype specified")
    end
end

function format_chemistry_string(reactants::Array{Symbol}, products::Array{Symbol})
    #=
    Given arrays of the reactants and products in a reaction,
    this small function constructs a chemistry reaction string.
    =#
    return string(join(reactants, " + ")) * " --> " * string(join(products, " + "))
end

function format_scin(n)
    #=
    Return number n in tidy scientific notation with only one decimal
    =#
    return @sprintf "%.1E" n
end

function format_sec_or_min(t) 
    #=
    Takes a time value t in seconds and formats it to read as H,M,S.
    =#

    thehour = div(t, 3600.)
    themin = div(t%3600., 60.)
    thesec = round(t%3600. - (themin*60.), digits=1)

    return "$(thehour) hours, $(themin) minutes, $(thesec) seconds"
end

function string_to_latexstr(a; dollarsigns=true)
    #=
    Given some chemistry reaction string with things like "pl" and un-subscripted numbers, 
    this will format it as a latex string for easy plotting.
    =#
    if dollarsigns==true
        replaced = replace(a, "O1D"=>"O(^1D)", "Nup2D"=>"N(^2D)", "2pl"=>"_2^+", "3pl"=>"_3^+", "2"=>"_2", "3"=>"_3", "E"=>"e^-", 
                                      "J"=>"", "plp"=>"^+ +", "pl"=>"^+",  "p"=>"+", 
                                      "-->"=>"\\rightarrow", "to"=>" \\rightarrow ")
        
        returnme = latexstring("\$\\mathrm{" * replaced * "}\$" ) 
    else # This is if you're trying to interpolate a variable that contains a latex string into an existing latex string...
        replaced = replace(a, "1D"=>"^1D", "up2D"=>"^2D", "2"=>"_2", "3"=>"_3", "E"=>"e^-", "pl"=>"^+", "-->"=>"\\rightarrow")
        returnme = "\\mathrm{" * replaced * "}" 
    end

    return returnme
end




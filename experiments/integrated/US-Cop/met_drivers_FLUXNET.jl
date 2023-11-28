
global main_dir
global sites

function replace_with_longterm_mean!(tmp::Matrix, selected_columns::Vector{String}, actual_column_names::Vector{String},plot_flag::Bool)

    # Loop through each selected column
    for col_name in selected_columns
        # println(col_name)
        if col_name .== "GPP_DT_VUT_REF"
            col_name_QC = "NEE_VUT_REF_QC"
        elseif col_name .== "LE_CORR"
            col_name_QC = "LE_F_MDS_QC"

        else
            col_name_QC = col_name * "_QC"
        end
        col_idx = findfirst(==(col_name), actual_column_names)
        #if col_idx is empty
        if col_idx === nothing
        else
            col_idx_QC = findfirst(==(col_name_QC), actual_column_names)
            if col_name .== "SW_OUT"
            else
                nonzero_indices = findall(x -> x != 0, tmp[2:end, col_idx_QC]) .+ 1
                tmp[nonzero_indices, col_idx] .= missing
            end

            # Extract the column data, skipping the name in the first row
            col_data = tmp[2:end, col_idx]
            formatted_dates = DateTime.(string.(tmp[2:end, 1]), "yyyymmddHHMM")
            if (col_name .== "GPP_DT_VUT_REF") || (col_name .== "SWC_F_MDS_1") 
                indices_below_900 = findall(x -> ismissing(x) ? false : x < 0, col_data) #for GPP 
            else
                indices_below_900 = findall(x -> ismissing(x) ? false : x < -900, col_data)
            end
            col_data[indices_below_900].=missing
            valid_values = filter(x -> !(x === missing || isnan(x) || x < -900), col_data)

            if isempty(valid_values)
                if col_name .== "LW_IN_F"
                    col_name = "LW_IN_JSB"
                    col_name_QC = col_name * "_QC"
                    col_idx = findfirst(==(col_name), actual_column_names)
                    if col_idx === nothing
                        continue
                    else
                        println("$col_name is being used")
                        col_idx_QC = findfirst(==(col_name_QC), actual_column_names)
                        nonzero_indices = findall(x -> x != 0, tmp[2:end, col_idx_QC]) .+ 1
                        tmp[nonzero_indices, col_idx] .= missing
                        # Extract the column data, skipping the name in the first row
                        col_data = tmp[2:end, col_idx]
                        valid_values = filter(x -> !(x === missing || isnan(x) || x < -900), col_data)
                        if isempty(valid_values)
                            @error "considering all record data $col_name is not available"
                            continue
                        end
                        if plot_flag
                            global savedir_input
                            filename = joinpath(savedir_input, string(tmp[1, col_idx]) * ".png")
                            plot(tmp[2:end, col_idx], title=tmp[1, col_idx], dpi=400)
                            savefig(filename)
                            end
                    end
                end
                if col_name .== "SWC_F_MDS_1"
                    col_name = "SWC_F_MDS_2"
                    col_name_QC = col_name * "_QC"
                    col_idx = findfirst(==(col_name), actual_column_names)
                    #col_idx = findfirst(x -> occursin(col_name, x), actual_column_names)
                    if col_idx === nothing
                        continue
                    else
                        println("$col_name is being used")
                        col_idx_QC = findfirst(==(col_name_QC), actual_column_names)
                        nonzero_indices = findall(x -> x != 0, tmp[2:end, col_idx_QC]) .+ 1
                        indices_below_900 = findall(x -> ismissing(x) ? false : x < 0, col_data) #for SWC negative values do not bear any meaning
                        tmp[indices_below_900, col_idx] .= missing
                        tmp[nonzero_indices, col_idx] .= missing
                        # Extract the column data, skipping the name in the first row
                        col_data = tmp[2:end, col_idx]
                        formatted_dates = DateTime.(string.(tmp[2:end, 1]), "yyyymmddHHMM")
                        # Compute long-term mean excluding NaN, missing, or values less than -900
                        # Create a DataFrame
                        
                        valid_values = filter(x -> !(x === missing || isnan(x) || x < -900), col_data)
                        if isempty(valid_values)
                            @error "considering all record data $col_name is not available"
                            continue
                        end
                        if plot_flag
                            global savedir_input
                            filename = joinpath(savedir_input, string(tmp[1, col_idx]) * ".png")
                            plot(tmp[2:end, col_idx], title=tmp[1, col_idx], dpi=400)
                            savefig(filename)
                            end
                    end
                end
                if col_name .== "TS_F_MDS_1"
                    col_name = "TS_F_MDS_2"
                  
                    col_name_QC = col_name * "_QC"
                    col_idx = findfirst(==(col_name), actual_column_names)
                    #col_idx = findfirst(x -> occursin(col_name, x), actual_column_names)
                    if col_idx === nothing
                        continue
                    else
                        println("$col_name is being used")
                        col_idx_QC = findfirst(==(col_name_QC), actual_column_names)
                        nonzero_indices = findall(x -> x != 0, tmp[2:end, col_idx_QC]) .+ 1
                        tmp[nonzero_indices, col_idx] .= missing
                        # Extract the column data, skipping the name in the first row
                        col_data = tmp[2:end, col_idx]
                        valid_values = filter(x -> !(x === missing || isnan(x) || x < -900), col_data)
                        if isempty(valid_values)
                            @error "considering all record data $col_name is not available"
                            continue
                        end
                        if plot_flag
                            global savedir_input
                            filename = joinpath(savedir_input, string(tmp[1, col_idx]) * ".png")
                            plot(tmp[2:end, col_idx], title=tmp[1, col_idx], dpi=400)
                            savefig(filename)
                            end
                    end
                else
                    continue
                end
            end
            formatted_dates = DateTime.(string.(tmp[2:end, 1]), "yyyymmddHHMM")
            # Compute long-term mean excluding NaN, missing, or values less than -900
            # Create a DataFrame
            df = DataFrame(date=formatted_dates, val=col_data)
            df[!, :week] = Dates.week.(df.date)
           # df[!, :day] = Day.(df.date.-Hour.(df.date).-df.date[1])
            df[!, :month] = month.(df[!, :date])
            # Filter the DataFrame in-place based on col_data values
            df_filtered = filter(row -> !(row.val === missing || isnan(row.val) || row.val < -900), df)
            grouped = groupby(df_filtered, :week)
            weekly_means = combine(grouped, :val => (x -> mean(skipmissing(x))) => :weekly_mean)
           # grouped = groupby(df_filtered, :day)
            #daily_mean = combine(grouped, :val => (x -> mean(skipmissing(x))) => :daily_mean)
            grouped = groupby(df_filtered, :month)
            monthly_mean = combine(grouped, :val => (x -> mean(skipmissing(x))) => :monthly_mean)
            # Merge the computed means back to the original DataFrame
            #df = leftjoin(df, daily_mean, on=:day)
            df = leftjoin(df, weekly_means, on=:week)
            df = leftjoin(df, monthly_mean, on=:month)
            #indices that are missing or nan values of -9999 are set to weekly mean
            indices = findall(x -> x === missing || isnan(x) || x < -900, df.val)
            all_indices = 1:length(col_data)
            # Find indices not in 'indices'
           # not_indices = setdiff(all_indices, indices)
            col_data[indices] = df.weekly_mean[indices]
            # the data that are present are replaced by daily mean
           # col_data[not_indices] = df.daily_mean[not_indices]
            #check again if there are any missing data left that was not filled with weekly data
            indices = findall(x -> x === missing || isnan(x) || x < -900, col_data)
            if .!isempty(indices)
                col_data[indices] = df.monthly_mean[indices]
                indices = findall(x -> x === missing || isnan(x) || x < -900, col_data)
                if .!isempty(indices)
                println("still some data is missing")
                end
            end
            # Put the updated data back into df
            tmp[2:end, col_idx] = col_data
            if plot_flag
            global savedir_input
            filename = joinpath(savedir_input, string(tmp[1, col_idx]) * ".png")
            plot(tmp[2:end, col_idx], title=tmp[1, col_idx], dpi=400)
            savefig(filename)
            end
        end
    end
end
function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    missing_indices = ismissing.(field)
    if sum(good_indices) .== 0# no good indices, then we use the data as it is
        field = vcat(field)[:]

        fill_value = mean(skipmissing(field))
        if sum(missing_indices) .== 0
        else
            field[missing_indices] .= fill_value
        end
    else
        fill_value = mean(skipmissing(field[good_indices]))
        field[.~good_indices] .= fill_value
        if sum(missing_indices) .== 0
        else
            field[missing_indices] .= fill_value
        end
        field = vcat(field)[:]
    end
    return field
end

function find_site_name_dirs(root_directory, site_name)
    directories_with_site_name = String[]

    for (root, dirs, files) in walkdir(root_directory)
        println(dirs)
        for dir in dirs

            if occursin(lowercase(site_name), lowercase.(dir))
                push!(directories_with_site_name, joinpath(root_directory, dir))
            end
        end
    end

    return directories_with_site_name
end
function find_csv_file_with_name(directory)
    for file in readdir(directory)
        path = joinpath(directory, file)
        if isfile(path) && occursin("FULLSET_HH", file) && endswith(file, ".csv")
            return path
        elseif isdir(path)
            result = find_csv_file_with_name(path)
            if result != nothing
                return result
            end
        end
    end
    return nothing
end
global site_name
global start_year
global end_year
global lat
global long
global spinup
#read fluxnet dataset and LAI from modis, if LAI is absent remove the corresponding dates from  FLUXNET sites
root_data_files = joinpath(metadata_dir, "FLUXNETSites/FLX_$site_name")
csv_file_path = find_csv_file_with_name(root_data_files)

driver_data_1 = readdlm(csv_file_path, ',')
column_names = string.(driver_data_1[1, :])
selected_columns = ["TA_F", "VPD_F", "PA_F", "P_F", "WS_F", "LW_IN_F", "SW_IN_F", "CO2_F_MDS", "SWC_F_MDS_1", "G_F_MDS", "H_F_MDS", "TS_F_MDS_1", "GPP_DT_VUT_REF", "LE_CORR", "H_F_MDS", "LW_OUT", "SW_OUT"]
#Stefan-Boltzmann 
replace_with_longterm_mean!(driver_data_1, selected_columns,column_names, plot_input_var)#all missing and below -900 datasets are replaced by long-term weekly mean
dates = driver_data_1[2:end, column_names.=="TIMESTAMP_START"]; # time yyyymmddhhmm
time_datadriver=DateTime.(string.(dates), "yyyymmddHHMM")
# LAI data
lai_dataset_path = joinpath(metadata_dir, "LAI");
files = readdir(lai_dataset_path)
lai_raw_data = joinpath(lai_dataset_path, files[occursin.(site_name, files)][1]);
global LAI_data = CSV.File(lai_raw_data) |> DataFrame
tmp1 = string.(LAI_data[:, "TIMESTAMP"])
tmptime = map(x -> DateTime(x, "yyyymmddHH"), tmp1)
#start
start_time_datadriver = minimum(time_datadriver)
end_time_datadriver = maximum(time_datadriver)
#end
start_tmptime = minimum(tmptime)
end_tmptime = maximum(tmptime)

# Determine the overlap range
overlap_start = max(start_time_datadriver, start_tmptime)
overlap_end = min(end_time_datadriver, end_tmptime)
filtered_time_datadriver = filter(x -> x >= overlap_start && x <= overlap_end, time_datadriver)
#filter 
(m,n)=size(driver_data_1)
matched_time=Vector{Any}(undef, m)
matched_time[1]=1
matched_time[2:end] = [time in filtered_time_datadriver for time in time_datadriver]
m=sum(matched_time)
driver_data=Array{Any}(undef, m, n)
driver_data=driver_data_1[matched_time .== 1, :]
date_strings = vec([string(date) for date in dates])
shortened_date_strings = [date[1:4] for date in date_strings]
year = parse.(FT, shortened_date_strings)
unique_year = unique(year)
first_year = unique_year[start_year]
last_year = unique_year[end_year]
filtered_indices = findall(year .>= first_year .&& year .<= last_year) .+ 1;
dates = driver_data[filtered_indices, column_names.=="TIMESTAMP_START"]; #
filtered_dates_date = DateTime.(string.(dates), "yyyymmddHHMM")

first_date = Dates.format(filtered_dates_date[1], "yyyy-mm-dd-HH")#"2004-01-01-06"

if spinup .== 1
    global nspinup = 1#number of years spinup occurs
    spinup_indices = findall(year .>= year[1] .&& year .<= year[1+nspinup]) .+ 1
    spinup_indices_dates = driver_data[spinup_indices, column_names.=="TIMESTAMP_START"] #
    spinup_dates = DateTime.(string.(spinup_indices_dates), "yyyymmddHHMM")
    driver_data = vcat(driver_data[1:1, :], driver_data[spinup_indices, :], driver_data[filtered_indices, :])
    global N_spinup_days = FT(365)
    global t_spinup = FT(365 * 24 * 60 * 60)#day to seconds

else
    global N_spinup_days = 0
    filtered_indices = [1; filtered_indices]
    driver_data = driver_data[filtered_indices, :]#only selected range of data is used
end

#replace here with columname final
TA = driver_data[2:end, column_names.=="TA_F"] .+ 273.15; # convert C to K
VPD = driver_data[2:end, column_names.=="VPD_F"] .* 100; # convert hPa to Pa
PA = driver_data[2:end, column_names.=="PA_F"] .* 1000; # convert kPa to Pa
P = driver_data[2:end, column_names.=="P_F"] ./ (1000 * 1800); # convert mm/HH to m/s
WS = driver_data[2:end, column_names.=="WS_F"]; # already m/s
LW_IN = driver_data[2:end, column_names.=="LW_IN_JSB"]
LW_IN_QC = driver_data[2:end, column_names.=="LW_IN_JSB_QC"]
replace_missing_with_mean!(LW_IN, LW_IN_QC)
SW_IN = driver_data[2:end, column_names.=="SW_IN_F"]
CO2_F = driver_data[2:end, column_names.=="CO2_F_MDS_QC"]
CO2 = driver_data[2:end, column_names.=="CO2_F_MDS"] .* 1e-6; # convert \mumol to mol
replace_missing_with_mean!(CO2, CO2_F)
SWC_F = driver_data[2:end, column_names.=="SWC_F_MDS_2_QC"] # Most likely 5cm depth
global SWC = driver_data[2:end, column_names.=="SWC_F_MDS_2"] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC, SWC_F)
TS_F = driver_data[2:end, column_names.=="TS_F_MDS_2_QC"] # Most likely 5cm depth
TS = driver_data[2:end, column_names.=="TS_F_MDS_2"] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS, TS_F)
global GPP = driver_data[2:end, column_names.=="GPP_DT_VUT_REF"] .* 1e-6 # to convert from micromol to mol.
LE = driver_data[2:end, column_names.=="LE_CORR"]
H_CORR = driver_data[2:end, column_names.=="H_CORR"]
H = driver_data[2:end, column_names.=="H_F_MDS"]
H_F = driver_data[2:end, column_names.=="H_F_MDS_QC"]
replace_missing_with_mean!(H, H_F)
G = driver_data[2:end, column_names.=="G_F_MDS"]
G_F = driver_data[2:end, column_names.=="G_F_MDS_QC"]
replace_missing_with_mean!(G, G_F)

LW_OUT = driver_data[2:end, column_names.=="LW_OUT"]# This has missing data
SW_OUT = driver_data[2:end, column_names.=="SW_OUT"]# This has missing data
#replace_missing_with_mean_by_value!(SW_OUT)
#replace_missing_with_mean_by_value!(LW_OUT)

global LOCAL_DATETIME = DateTime.(string.(driver_data[2:end, 1]), "yyyymmddHHMM")
 UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(6)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds
if spinup .== 0
    global t_spinup = DATA_DT
end
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )
e = @. esat - VPD
q = @. 0.622 * e ./ (PA - 0.378 * e)

#Make a bunch of splines
seconds = FT.(0:DATA_DT:((length(UTC_DATETIME)-1)*DATA_DT));
p_spline = Spline1D(seconds, -P[:]) # m/s
atmos_q = Spline1D(seconds, q[:])
atmos_T = Spline1D(seconds, TA[:])
atmos_p = Spline1D(seconds, PA[:])
atmos_co2 = Spline1D(seconds, CO2[:])
atmos_u = Spline1D(seconds, WS[:])
LW_IN_spline = Spline1D(seconds, LW_IN[:])
SW_IN_spline = Spline1D(seconds, SW_IN[:])
atmos_h = FT(32)
precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = eltype(t)(0) # this is likely not correct


# Construct the drivers
atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    UTC_DATETIME[1],
    atmos_h;
    c_co2=atmos_co2
)
#lat = FT(38.7441) # degree
#long = FT(-92.2000) # degree
datetime_obj = DateTime(first_date, "yyyy-mm-dd-HH")

# Add 8 hours
new_datetime = datetime_obj + Hour(8)#it should be UTC time, so for USA sites it is +8 hours

function zenith_angle(
    t::FT,
    orbital_data,
    ref_time;
    latitude=lat,
    longitude=long,
    insol_params=earth_param_set.insol_params
) where {FT}
    # This should be time in UTC
    dt = ref_time + Dates.Second(round(t))
    FT(
        instantaneous_zenith_angle(
            dt,
            orbital_data,
            longitude,
            latitude,
            insol_params,
        )[1],
    )
end

radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline,
    UTC_DATETIME[1];
    θs=zenith_angle,
    orbital_data=Insolation.OrbitalData()
)



tmp = [LAI_data[:, :Vcmax25_leaf_v1]; LAI_data[:, :Vcmax25_leaf_v2]];
nan_vec = replace(tmp, "NA" => "NaN") #global_clumping_index_uncompressed.tif
parsed_values = parse.(Float64, nan_vec)
parsed_values = [isnan(val) ? missing : val for val in parsed_values]
global parse_values
nan_vec = replace(LAI_data[:, "LAI"], "NA" => "NaN")
LAI_data = DataFrame([tmptime, parse.(Float64, nan_vec)], [:time_l, :LAI])
#take monthly average of LAI across all years and use it to gapfill the data
LAI_data[!, :month] = month.(LAI_data.time_l)
monthly_means = combine(groupby(LAI_data, :month)) do sdf
    mean_value = mean(filter(x -> !isnan(x), sdf.LAI))
    DataFrame(mean_LAI=mean_value)
end
#join back mean_LAI on the LAI data
LAI_data = leftjoin(LAI_data, monthly_means, on=:month)
# Replace NaN in LAI with the corresponding monthly mean
LAI_data[!, :LAI] = ifelse.(isnan.(LAI_data[!, :LAI]), LAI_data[!, :mean_LAI], LAI_data[!, :LAI])
mask = [t in filtered_dates_date for t in LAI_data.time_l]
LAI_filtered = LAI_data[mask, "LAI"]

if spinup .== 1
    spinup_indices_dates = DateTime.(string.(spinup_indices_dates), "yyyymmddHHMM")
    mask2 = [t in spinup_dates for t in LAI_data.time_l]
    spinup_LAI = LAI_data[mask2, "LAI"]
    LAI_repeated = vcat(repeat(spinup_LAI, inner= Int64(length(seconds)./365)), repeat(LAI_filtered, inner=48))
else
    LAI_repeated = repeat(LAI_filtered, inner= Int64(length(seconds)./365))
end
plot(LAI_filtered, title="LAI", dpi=400)
filename = joinpath(savedir_input,  "LAI_org.png")
savefig(filename)
plot(round.(seconds./(3600*24),digits=2),LAI_repeated, title="LAI repeated", dpi=400)
filename = joinpath(savedir_input,  "LAI_repeated.png")
savefig(filename)
# This has the same timestamp as the driver data, so it's ok to use the time column from that file here
LAIspline = Spline1D(seconds, LAI_repeated[:])
LAIfunction = (t) -> eltype(t)(LAIspline(t))

import cdsapi

dataset = "reanalysis-era5-pressure-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "geopotential",
        "u_component_of_wind",
        "v_component_of_wind",
	"temperature",
    ],
    "year": ["2026"],
    "month": ["03"],
    "day": ["05"],
    "time": ["12:00"],
    "pressure_level": ["500", "700", "900"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [90, -60, 0, 90] # up left down right
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


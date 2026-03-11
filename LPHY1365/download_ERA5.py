import cdsapi

dataset = "reanalysis-era5-pressure-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "geopotential",
        "u_component_of_wind",
        "v_component_of_wind"
    ],
    "year": ["2026"],
    "month": ["03"],
    "day": ["05"],
    "time": ["12:00"],
    "pressure_level": ["500", "1000"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [80, -60, 20, 70]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


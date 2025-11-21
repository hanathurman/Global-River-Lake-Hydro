#import collections
#try:
#    collections.MutableMapping = collections.abc.MutableMapping
#except AttributeError:
#    pass
#import zarr

import numpy as np
import pandas as pd
import s3fs
import xarray
import datetime


def pull_tributary(reach_id, start_date):
    start_date = datetime.datetime.strptime(start_date, "%m-%d-%Y").date()
    end_date = datetime.date.today()
    sim_begin = datetime.date(1940, 1, 1)
    first_index = start_date - sim_begin
    second_index = end_date - sim_begin
    bucket_uri = "s3://geoglows-v2/retrospective/daily.zarr"
    region_name = "us-west-2"
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name=region_name))
    s3store = s3fs.S3Map(root=bucket_uri, s3=s3, check=False)
    ds = xarray.open_zarr(s3store)
    df = (
        ds["Q"]
        .sel(river_id=reach_id, time=slice(start_date, end_date))
        .to_dataframe()
    )

    df = df.reset_index().set_index("time").pivot(columns="river_id", values="Q")
    return df

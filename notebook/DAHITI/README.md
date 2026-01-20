To use the notebook [_download_API_v2.ipynb_](./download_API_v2.ipynb) you need to files in the same directory as the notebook:

1. A text file ([_api_key.txt_](./api_key.txt)) containing your personal API key. Once you're logged in, you can request your API key in this [link](https://dahiti.dgfi.tum.de/en/frequently-asked-questions/api-key/).
2. A YML file ([config.yml](./config.yml)) with the arguments used to filter the targets and variables of interest. You can select targets based on country and type of point:
    * `api_url`: the URL to which the request will be submitted. By default: https://dahiti.dgfi.tum.de/api/v2/
    * `type`: kind of points of interest (lake, reservoir, river or wetland)
    * `country`: a list of country of interst
    * `variables`: a list of variables of interest (water-level, surface-area, volume-variation, hypsometry). By default, only water-level.
    * `ouput_path`: rood directory where the data will be saved. The code will create separate subfolders for the _targets_ and the _time_series_.
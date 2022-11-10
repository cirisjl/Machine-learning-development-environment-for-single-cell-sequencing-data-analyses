import requests
import os
import os.path as osp
import aiofiles
import aiohttp
import errno
# from auth.credentials import PATH

def download_file_from_google_drive(destination):
    #URL = "https://docs.google.com/uc?export=download"
    URL="https://ftp.cdc.gov/pub/health_statistics/nchs/Publications/ICD10CM/2017/icd10cm_codes_2017.txt"

    session = requests.Session()

    # response = session.get(URL, params = { 'id' : id }, stream = True)
    response=session.get(URL,stream=True)
    print(str(response))
    token = get_confirm_token(response)
    print(token)
    if token:
        # params = { 'id' : id, 'confirm' : token }
        params={'confirm':token}
        response = session.get(URL, params = params, stream = True)
    else:
        params={}
    print(params)
    save_response_content(response, destination)    

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
def makedirs(path):
    try:
        os.makedirs(osp.expanduser(osp.normpath(path)))
    except OSError as e:
        if e.errno != errno.EEXIST and osp.isdir(path):
            raise e

async def async_download(url, folder, log=True):
    r"""Downloads the content of an URL to a specific folder.
    Args:
        url (string): The url.
        folder (string): The folder.
        log (bool, optional): If :obj:`False`, will not print anything to the
            console. (default: :obj:`True`)
    """
    
    filename = url.rpartition('/')[2]
    path = osp.join(folder, filename)
    
    if osp.exists(path) and osp.getsize(path) > 0:  # pragma: no cover
        if log:
            print('Using exist file', filename)
        return path

    if log:
        print('Downloading', url)

    makedirs(folder)
    print(makedirs)
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:

                chunk_size = 1024*1024

                downloaded_size = 0
                
                async with aiofiles.open(path, 'wb') as f:
                        
                    async for data in response.content.iter_chunked(chunk_size):
                        downloaded_size += len(data)
                        await f.write(data)
                    
                    return path
    except:
        if os.path.exists(path):
             os.remove(path)
        raise RuntimeError('Stopped downloading due to interruption.')

if __name__ == "__main__":
    # file_id = '1TCy9OoCRkNWzo22fGbycG-j9z7gphfGE'
    path = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/"
    destination = path + "tmp/hmi.txt"
    destination2=path+"tmp/"
    download_file_from_google_drive(destination)
    # await async_download(url='https://drive.google.com/file/d/1TCy9OoCRkNWzo22fGbycG-j9z7gphfGE/view?usp=sharing',folder=destination2,log=True)
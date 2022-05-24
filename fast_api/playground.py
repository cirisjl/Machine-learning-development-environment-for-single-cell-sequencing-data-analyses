import aiohttp
import asyncio
from tqdm import tqdm
import aiofiles
import os.path as osp

GBFACTOR = float(1 << 30)

async def main():

    url = 'http://snap.stanford.edu/ogb/data/graphproppred/csv_mol_download/bace.zip'
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:

            chunk_size = 1024*1024

            downloaded_size = 0
            
            folder = "/Users/jsaied/Library/Mobile Documents/com~apple~CloudDocs/University/GS/obj/fast-api/tmp"
            
            filename = url.rpartition('/')[2]
            path = osp.join(folder, filename)
            
            async with aiofiles.open(path, 'wb') as f:
                    
                async for data in response.content.iter_chunked(chunk_size):
                    downloaded_size += len(data)
                    await f.write(data)


loop = asyncio.get_event_loop()
loop.run_until_complete(main())
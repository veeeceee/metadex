import aiofiles
import aiohttp
import asyncio
import async_timeout
import os
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

def create_path(path):
    """
    Summary: creates a path, with error checking

    Args:
        path (str): name of folder/directory to be created

    Returns:
        None, creates folder
    """

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

async def download_coroutine(session, url, group, sample, source, evalue, identity, length):
    """
    Summary:

    Args:

    Returns:
    """

    async with session.get(url, params = {'source': source,'evalue': evalue, 'identity': identity, 'length': length, 'type': 'all'}, timeout=None) as response:
        filename = '%(group)s_%(sample)s.tsv' % \
            {"group": group, "sample": sample}
        print('reading '+ str(filename))

        async with aiofiles.open(filename, 'wb') as fd:
            print('writing ' + str(filename))
            while True:
                chunk = await response.content.read(1024)
                if not chunk:
                    break
                await fd.write(chunk)
        return await response.release()


async def download_geneonly_coroutine(session, url, study, idList, source, evalue, identity, length):
    """
    Summary:

    Args:

    Returns:
    """
    params_ids = {'id': mgID for mgID in idList}
    params_other =  {'source': source,'evalue': evalue, 'identity': identity, 'length': length}
    params = params_ids.copy()
    params.update(params_other)
    geneOnlyURL = 'http://api.metagenomics.anl.gov/matrix/function?id=' + '&id='.join(idList)+ '&source=' + str(source)+ '&evalue='  + str(evalue) + '&identity=' + str(identity) + '&length=' + str(length)

    async with session.get(url, params = params, timeout=None) as response:
        filename = '%(study)s_function_counts.biom' % \
            {"study": study}

        async with aiofiles.open(filename, 'wb') as fd:
            while True:
                chunk = await response.content.read(128)
                if not chunk:
                    break
                await fd.write(chunk)
        return await response.release()





async def download_all(urls, geneOnlyUrl, study, metagenomeGroupDF, source, evalue, identity, length):
    """Launch requests for all web pages."""
    tasks = []
    idList = metagenomeGroupDF['metagenome ID']
    groupList = metagenomeGroupDF['group']
    sampleList = metagenomeGroupDF['sample']
    create_path(str(study))
    os.chdir(str(study))
    async with aiohttp.ClientSession() as session:
        for i in range(len(urls)):
            task = asyncio.ensure_future(download_coroutine(session, urls[i], groupList[i], sampleList[i], source, 5, 60, 15))
            tasks.append(task) # create list of tasks
        tasks.append(asyncio.ensure_future(download_geneonly_coroutine(session, geneOnlyUrl, study, idList, source, 5, 60, 15)))
        _ = await asyncio.gather(*tasks) # gather task responses


def get_all_async(study, metagenomeGroupDict, source, evalue=5, identity=60, length=15):
    """
    Summary: gets all metagenomic annotations for a study from the MG-RAST servers and stores as CSV

    Args:
        study (str)
        metagenomeGroupDict (dict)
        source (str)
        evalue (int)
        identity (int)
        length (int)

    Returns:
        None, saves all annotations in folder named study
    """
    metagenomeGroupDF = pd.DataFrame({'metagenome ID' : list(metagenomeGroupDict.keys()),'group' : list(metagenomeGroupDict.values())})
    metagenomeGroupDF['sample'] = metagenomeGroupDF.groupby('group').cumcount()+1
    urls = ["http://api.metagenomics.anl.gov/annotation/similarity/" + str(id) for id in metagenomeGroupDF['metagenome ID']]
    geneOnlyUrl = 'http://api.metagenomics.anl.gov/matrix/function'
    loop = asyncio.get_event_loop() # event loop
    future = asyncio.ensure_future(download_all(urls, geneOnlyUrl, study, metagenomeGroupDF, source, evalue, identity, length)) # tasks to do
    loop.run_until_complete(future) # loop until done
    os.chdir('..')


if __name__ == '__main__':
    get_all_async('lagoon study', {'mgm4739968.3':'nador_lagoon', 'mgm4739969.3':'oualidia_lagoon', 'mgm4739970.3': 'oualidia_lagoon', 'mgm4739971.3':'nador_lagoon'}, 'RefSeq', evalue=5, identity=60, length=15)

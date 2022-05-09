import re
import requests as req
import os

def get_html(url):
    r = req.get(url, params=None)
    text = r.text
    return text

def get_catalog_urls(text):
    regex = r'(?<=<a\shref=\')(.*?)(?=\')'
    filter_regex = r'catalog.*.nc'
    URLS = re.findall(regex, text)
    filtered = []
    for i in URLS:
        if re.match(filter_regex, i):
            filtered.append(i)
    return filtered

def urls_to_data(url):
    html_body = get_html(url)
    filtered = get_catalog_urls(html_body)
    mean_regex = r'catalog.*mean.*.nc'
    clean_regex = r'catalog\.html\?dataset\='
    files_regex = r'_files'
    rep_str  = 'https://thredds.met.no/thredds/dodsC/fou-hi/'
    no_mean = []

    #filtering out mean data
    for i in filtered:
        if not re.match(mean_regex, i):
            no_mean.append(i)
    
    #changing name to thredds url and removing "_files"
    for i, j in enumerate(no_mean):
        tmp = re.sub(clean_regex, rep_str, j)
        no_mean[i] = re.sub(files_regex, '', tmp)
    
    return no_mean

def save_to_file(list):
    out = open('thredds_urls.txt','w')
    for i in list:
        out.write(f'{i}\n')


if __name__ == '__main__':
    urls = urls_to_data('https://thredds.met.no/thredds/catalog/fou-hi/barents_eps_eps/catalog.html')
    save_to_file(urls)
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps/barents_epsmean_20220421T18Z.nc.html
#catalog.html?dataset=barents_eps_eps_files/barents_eps_20220419T18Z.nc
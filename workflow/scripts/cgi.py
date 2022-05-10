# -*- coding: utf-8 -*-

email = "jouni.kujala@uef.fi"
token = "7d7d1252c6136a62d0db"


import requests
import time

def run_cgi(email, token):
    
    url = "https://www.cancergenomeinterpreter.org/api/v1/"

    # Define used header.
    headers = {"Authorization": "{} {}".format(email, token)}
    
    # Define analysis parameters.
    payload = {"cancer_type": "BRCA", "reference": "hg38"}
    
    # Define input files.
    files = {"mutations": open("Mer3_filtered_vep.vcf", "rb")}

    # Start analysis.
    request = requests.post(url, headers=headers, files=files, data=payload)
    job_id = request.json()
    
    status = "Running"
    payload={"action":"logs"}
    
    while status != "Done":
        time.sleep(10)
        print("Running...")
        
        status = requests.get(url + job_id, headers=headers, params=payload)
        status = status.json()["status"]

    # Download results.
    payload={"action":"download"}
    request = requests.get(url + job_id, headers=headers, params=payload)
    
    with open('file.txt', 'wb') as fd:
        fd.write(request._content)

    # Delete analysis.
    print("Removing analysis...")
    requests.delete(url + job_id, headers=headers)

run_cgi(email, token)

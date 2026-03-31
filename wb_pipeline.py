!pip install biopython
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re

Entrez.email = "your_email@example.com"


# ----------------------------
# PubMed 검색
# ----------------------------
def search_pubmed(protein, max_results=20):

    query = f"{protein} western blot"

    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_results
    )

    record = Entrez.read(handle)

    return record["IdList"]


# ----------------------------
# PMID -> PMCID
# ----------------------------
def get_pmcid(pmid):

    url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={pmid}&format=json"

    r = requests.get(url).json()

    try:
        return r["records"][0]["pmcid"]
    except:
        return None


# ----------------------------
# PMC XML 다운로드
# ----------------------------
def download_xml(pmcid):

    url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/?format=xml"

    r = requests.get(url)

    if r.status_code == 200:
        return r.text

    return None


# ----------------------------
# figure caption 추출
# ----------------------------
def extract_figure_captions(xml_text):

    soup = BeautifulSoup(xml_text, "lxml")

    captions = []

    for fig in soup.find_all("fig"):

        cap = fig.find("caption")

        if cap:
            captions.append(cap.get_text(" ", strip=True))

    return captions


# ----------------------------
# 본문 문장 추출
# ----------------------------
def extract_body_sentences(xml_text):

    soup = BeautifulSoup(xml_text, "lxml")

    sentences = []

    for p in soup.find_all("p"):

        text = p.get_text(" ", strip=True)

        split = re.split(r'\.|\n', text)

        sentences.extend(split)

    return sentences


# ----------------------------
# western blot 필터
# ----------------------------
def filter_wb_sentences(sentences):

    wb_sentences = []

    for s in sentences:

        s_lower = s.lower()

        if any(k in s_lower for k in [
            "western blot",
            "immunoblot",
            "wb analysis"
        ]):

            wb_sentences.append(s.strip())

    return wb_sentences


# ----------------------------
# effect detection
# ----------------------------
def detect_effect(sentence):

    s = sentence.lower()

    if any(w in s for w in [
        "increase",
        "upregulate",
        "elevated",
        "enhanced"
    ]):
        return "increase"

    if any(w in s for w in [
        "decrease",
        "reduce",
        "downregulate",
        "suppressed"
    ]):
        return "decrease"

    return "unknown"


# ----------------------------
# phospho detection
# ----------------------------
def detect_phospho(sentence):

    s = sentence.lower()

    if any(k in s for k in [
        "phosphorylation",
        "phospho",
        "p-"
    ]):
        return "phospho"

    return "total"


# ----------------------------
# gene knockdown detection
# ----------------------------
def detect_knockdown(sentence):

    s = sentence.lower()

    if any(k in s for k in [
        "knockdown",
        "sirna",
        "shrna",
        "silencing"
    ]):
        return True

    return False


# ----------------------------
# drug detection (사용자 입력 반영)
# ----------------------------
def detect_drug(sentence, user_drug=None):

    s = sentence.lower()

    # 사용자 입력 drug 우선
    if user_drug:
        if user_drug.lower() in s:
            return user_drug
        else:
            return None

    # 기본 drug 리스트 (optional)
    drug_list = [
        "ly294002",
        "wortmannin",
        "rapamycin",
        "torin",
        "mk2206"
    ]

    for d in drug_list:
        if d in s:
            return d

    return None


# ----------------------------
# 메인 파이프라인
# ----------------------------
def run_pipeline(protein, drug=None, max_results=20):

    pmids = search_pubmed(protein, max_results)

    results = []

    for pmid in pmids:

        pmcid = get_pmcid(pmid)

        if not pmcid:
            continue

        xml = download_xml(pmcid)

        if not xml:
            continue

        captions = extract_figure_captions(xml)
        body = extract_body_sentences(xml)

        sentences = captions + body

        wb_sentences = filter_wb_sentences(sentences)

        for s in wb_sentences:

            if protein.lower() not in s.lower():
                continue

            effect = detect_effect(s)
            phospho = detect_phospho(s)
            kd = detect_knockdown(s)
            drug_found = detect_drug(s, drug)

            results.append({
                "PMID": pmid,
                "protein": protein,
                "form": phospho,
                "drug": drug_found,
                "gene_KD": kd,
                "effect": effect,
                "sentence": s
            })

    return pd.DataFrame(results)


# ----------------------------
# 실행 (🔥 사용자 입력)
# ----------------------------

protein = input("Enter protein name: ")
drug = input("Enter drug name (optional, press Enter to skip): ")

if drug == "":
    drug = None

df = run_pipeline(protein, drug, max_results=30)

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

print(df.head())

# 파일명 자동 생성
filename = f"{protein}_western_blot_results.csv"
df.to_csv(filename, index=False)

print(f"\nSaved results to {filename}")

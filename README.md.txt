# Western Blot Text Mining Pipeline

## Overview
This project automatically extracts Western blot-related information from PubMed/PMC articles.

## Features
- PubMed search via Entrez API
- PMC XML parsing
- Western blot sentence detection
- Effect classification (increase/decrease)
- Phosphorylation detection
- Drug & gene knockdown detection

## Tech Stack
- Python
- BioPython (Entrez)
- BeautifulSoup
- Pandas

##  Example
```python
df = run_pipeline("AKT", max_results=30)
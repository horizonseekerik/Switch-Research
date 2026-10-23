import sys
sys.path.append(r"C:\Users\hp\.gemini\config\plugins\science\skills\literature_search_arxiv\scripts")
from search_arxiv import _CLIENT
import urllib.parse
import xml.etree.ElementTree as ET
import json

q = 'all:"phase change" AND all:"optical switch" AND all:silicon'
params = {'search_query': q, 'start': 0, 'max_results': 8, 'sortBy': 'relevance', 'sortOrder': 'descending'}
url = 'http://export.arxiv.org/api/query?' + urllib.parse.urlencode(params)
xml_data = _CLIENT.fetch_bytes(url)
root = ET.fromstring(xml_data)
ns = {'atom': 'http://www.w3.org/2005/Atom'}
papers = []
for entry in root.findall('atom:entry', ns):
    title = entry.find('atom:title', ns).text.strip().replace('\n', ' ')
    summary = entry.find('atom:summary', ns).text.strip().replace('\n', ' ')
    pdf_url = entry.find('atom:id', ns).text.strip()
    papers.append({'title': title, 'summary': summary, 'url': pdf_url})

with open('arxiv_broad_papers.json', 'w', encoding='utf-8') as f:
    json.dump(papers, f, indent=2)
print(f"Retrieved {len(papers)} papers.")

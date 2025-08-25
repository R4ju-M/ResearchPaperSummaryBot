import pandas as pd
from openai import OpenAI
from dotenv import load_dotenv
import os
import smtplib, ssl
from Bio import Entrez
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup, Comment
import requests

#Creating Open API Key to allow access to GPT
#Key stored in .env file
load_dotenv()
client = OpenAI(
    api_key = os.getenv("OPEN_API_KEY")
)

#Provide your own email for security measures for Entrez
Entrez.email = "dal367730@utdallas.edu"
referenceSynonyms = ["reference", "references", "citations", "abbreviations"]

def extract_keywords(user_input, max_attempts = 3):
    attempt = 0
    #Prompt defining the OpenAI bots role in extracting keywords from the user query
    final_prompt = f""""You are a PubMed search assistant. Your job is to extract keywords from the user's research interest
        and format them into a PubMed-compatible search query. Use Boolean operators (AND, OR) appropriately.
        Only return the search query string. Do not include explanations or extra text. Use OR operator to include very specific synonyms such as acryonyms
        and similar things. Do not use broad phrases as synonyms like gene editing unless specifically mentioned.
        You are on attempt {attempt}. If you're on attempt = 0, be ultraspecific to the keywords mentioned, but still use tightly related synonyms/acronyms to improve search results. As it increases,
        get very slightly more specific to return somewhat relevant articles, IF NEED BE. 
        """
    
    #giving it 3 tries to give accurate keywords, so while attemps less than 3, loop will run
    while attempt < max_attempts: 
        completion = client.chat.completions.create(
        model = "gpt-4o-mini",
        #defining the role of the LLM agent as well as what the user will prompt it
        messages = [
                {
                "role": "system",
                "content": final_prompt
                },
                {
                "role": "user",
                "content": f"User query: {user_input}\n\n."
                }
            ]
        )

        #extracting output from completion module of OpenAI API
        kw = completion.choices[0].message.content.strip()
        print(f"[Attempt {attempt+1}] Extracted Query: {kw}")

        #Will return a maximum of 10 articles that contain the keywords
        handle = Entrez.esearch(db="pubmed", term=kw, retmax = 10)
        record = Entrez.read(handle)
        handle.close()

        #If articles are found, we are done. If no articles were found from the keywords then go onto the next attempt.
        #The next attempt will make the keywords broader to heighten chances of getting matching articles
        id_list = record.get("IdList", [])
        if id_list:
            print(f"Found {len(id_list)} articles.\n")
            return kw
        else: 
            print("No results found. Retrying.\n")
            attempt += 1

    print("No results found for given query.")
    return None

#method formats the idList provided in XML format into a loopable array
def getIdList(keywords, max_results = 3):
    handle = Entrez.esearch(db = "pmc", term = keywords, retmax = max_results)
    records = Entrez.read(handle)
    handle.close()

    idList = records.get('IdList', [])

    if not idList:
        return []
    else:
        return idList

#searches for comments made in the XML file of the article. We want to search for any comments that mention no text-mining allowed
def iscomment(elem):
    return isinstance(elem, Comment)

def getArticleTextBody(pubMedId): 
    pub_id = pubMedId
    print(pub_id)
    fullTextBody = ""
    handle = Entrez.efetch(db = "pmc", id = pub_id, retmode = "xml", rettype = "full")
    #extracting all the XML content of the article
    xml_data = handle.read()
    #using BeautifulSoup API tool to help filter through the XML file output for what we want - the article text
    soup = BeautifulSoup(xml_data, "xml")

    #check whether or not the article allows text-mining
    comments = soup.find_all(string=iscomment)
    for comment in comments:
        if "does not allow downloading of the full text" in comment.lower():
            print("Article does not allow text mining, but here is the PMC ID: " + pubMedId)
            return None
    
    articleSectionsList = []
    
    #Extracting only the text of the article
    #First relevant part of the article is the abstract, so we will skip to the abstract section
    abstract_tag = soup.find('abstract')
    if abstract_tag:
        fullTextBody += "*** Abstract ***\n"
        for p in abstract_tag.find_all('p'):
            fullTextBody += p.get_text() + "\n\n"
    
    #find all section tags, and then parse out the text from each section in the body (e.g. Methodology, Results, Conclusion, etc.)
    #go all the way down through the article till we hit the references section. It is full of citations that are not relevant to generating a summarry of the article
    for sec in soup.find_all('sec'):
        sectionText = ""
        trueTitle = ""
        title_tag = sec.find('title')
        
        if title_tag is not None:
            trueTitle = title_tag.get_text()
        else:
            trueTitle = "Untitled Section"
        
        if trueTitle.lower() in referenceSynonyms:
            print(trueTitle)
            break
        
        sectionText += "***" + trueTitle + "***\n"
        for paragraph_tag in sec.find_all('p'):
            paragraphText = paragraph_tag.get_text()
            sectionText += paragraphText + "\n\n"

        articleSectionsList.append(sectionText)


    fullTextBody += "".join(articleSectionsList)

    #This code can be commented out if you want
    #Just a way to see the entirety of the text of the article in a text file called "storage"
    with open("storage.txt", "w") as file:
        file.write("")
    with open("storage.txt", "a") as file:
        file.write(fullTextBody)
        file.write("\n\n\n___________________________________")
    
    return fullTextBody

def getSummaryOfArticle(articleTextBody, userInput):
    #defining role of API tool to create summaries of the research article we provide it
    systemPrompt = f"""
    You are an expert biomedical research assistant. You will be given a user query and the full text of a biomedical research article.  
    Read the article carefully and summarize only information directly relevant to the query.  
    Ignore unrelated sections and do not invent information. If the article does not mention the query, state it contains no relevant information.  
    Keep the summary concise and focused in a single paragraph. 
    """.strip()


    completion = client.chat.completions.create(
        model="gpt-4o-mini", 
        messages = [
        {
            "role": "system",
            "content": systemPrompt
        },
        {
            "role": "user",
            "content": (
                f"User Query: {userInput}\n\n"
                f"Article Text:\n{articleTextBody}\n\n"
                "Please generate a single-paragraph summary that only includes information directly relevant to the query. Ignore unrelated content."
            )
        }
    ],

        #Temperature ranges from 0-2. The closer to 0, the more conservative and focused the results will be. Higher temperatures will be more random and creative
        #Summary limited to 300 words for brevity
        temperature=0.3,
        max_tokens=300
    )

    return completion.choices[0].message.content.strip()

def getArticleTitle(pubMedId):
    raw_records = Entrez.efetch(db="pmc", id = pubMedId, rettype = "medline", retmode = "text")
    read_records = raw_records.read()
    raw_records.close()

    for line in read_records.split("\n"):
        #extracting title out of the list of information in entry
        if line.startswith("TI  - "):
            title = line.replace("TI  - ", "").strip()
            return title
        
    return "No Title Found"

#Enter query: E.g. Please give me summaries of articles that touch on the use of AAV particles to deliver CRISPR-Cas13 systems
userInput = input("Please enter your research query: ")

#Clears output.txt every time program is run
with open("output.txt", "w") as file:
    file.write("")

keywords = extract_keywords(userInput)
idList = getIdList(keywords)

#Iterates through each ID, providing the title and the corresponding, relevant summary of the article
for id in idList:
    articleTitle = getArticleTitle(id)
    articleTextBody = getArticleTextBody(id)
    summaryOfArticle = getSummaryOfArticle(articleTextBody, userInput)
    with open("output.txt", "a") as file:
        file.write("Title: " + articleTitle + "\n\n")
        file.write(summaryOfArticle)
        file.write("\n________________________________________________________\n")




    

    


    












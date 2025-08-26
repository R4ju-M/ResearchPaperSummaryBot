# ResearchPaperSummaryBot
Program that generates revenant summaries of PMC articles based on the user query

# Installation

## 1. Clone the Repository
```
git clone https://github.com/R4ju-M/ResearchPaperSummaryBot.git
cd ResearchPaperSummaryBot
```
## 2. Open the folder in VS Code:
Open VS Code → File → Open Folder → select ResearchPaperSummaryBot.  
All your code and the requirements.txt will already be inside.  
  
## 3. Create a virtual environment
```
python -m venv venv        -- try typing python3 instead of python, if facing issue
source venv/bin/activate   -- On macOS/Linux
venv\Scripts\activate      -- On Windows
```
  
## 4. Install Dependencies
```
pip install -r requirements.txt
```

## 5. Create OpenAPI Key
Create a .env file in your folder location
Inside this file, type this exact line. Replace your_key with your OpenAPI variable key
```
OPEN_API_KEY = your_key
```

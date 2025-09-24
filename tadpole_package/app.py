import streamlit as st
import streamlit.components.v1 as components
import traceback
import matplotlib.pyplot as plt
import graphviz
import glob
import sbol3
import os
import RNA
# Import using absolute rout so it can be used as a package and launch the streamlit app
from tadpole import *
from scipy.spatial.distance import pdist, squareform
import numpy as np
import shutil
import subprocess
import time
from io import StringIO
from contextlib import redirect_stdout
# from weasyprint import HTML # Weasyprint is typically for PDF generation and might require specific server setup not always available in Streamlit Cloud
import tempfile
import datetime
from collections import defaultdict
import base64 
import pandas as pd # Added pandas import as it's used later for dataframes

st.set_page_config(page_title="Tadpole", layout="wide", page_icon="tadpole_package/images/logo.png") # Updated page title

# ---------------- CUSTOM STYLES ----------------
st.markdown("""
    <style>
            
    /* Global Styles */
    html, body, [class*="css"] {
        font-family: 'Helvetica Neue', sans-serif;
        background-color: #f8f8f8; /* Light background for the entire app */
        color: #333333; /* Default text color */
    }

    /* Primary Color Palette based on #e52920 */
   :root {
    --bg-light-red: #fdeeee;      /* Very light red/pink for backgrounds */
    --soft-red: #fbdad9;          /* Light red/pink for inactive elements, light card borders */
    --medium-light-red: #f7bcbb;  /* Slightly darker pink */
    --light-accent-red: #f28f8c;  /* Medium pink/light red for hover states */
    --primary-red: #ea534e;       /* Adjusted to #ea534e as per new palette, was #e52920 */
    --dark-accent-red: #e62e25;   /* Slightly darker primary red */
    --dark-red: #e41513;          /* Darker red for button hover/active */
    --very-dark-red: #be1818;     /* Even darker red */
    --deep-red: #9d1915;          /* Dark red/maroon for strong accents/borders */
    --text-darker-red: #53110e;   /* Very dark red for text/headers */
    
    --soft-gray: #e0e0e0;         /* Soft gray for general backgrounds/borders (kept as is) */
    --medium-gray: #cccccc;       /* Medium gray for borders/dividers (kept as is) */
    --dark-gray: #1a1a1a;         /* General dark text/header color (kept as is) */
    --text-color-light: #333333;  /* Dark text on light backgrounds (kept as is) */
    --text-color-dark: #ffffff;   /* White text on dark primary colors (kept as is) */
    --card-bg: #ffffff;           /* White background for cards/panels (kept as is) */
    --border-color: var(--soft-gray); /* Default border color (kept as is) */
    --shadow-color: rgba(0, 0, 0, 0.1); /* Subtle shadow (kept as is) */
}
    }

    /* Headers - Ensure good contrast */
    h1, h2, h3, h4, h5, h6 {
        color: var(--text-darker-red); /* Using the darkest red for headers */
        font-weight: 600;
    }
    h1 { font-size: 2.5em; color: var(--primary-red); margin-bottom: 0.5em; } /* Keep H1 primary red */
    h2 { font-size: 2em; margin-bottom: 0.5em; }
    h3 { font-size: 1.5em; margin-bottom: 0.5em; }
    

    /* Main Container Padding */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        padding-left: 3rem;
        padding-right: 3rem;
    }

    /* Tabs Styling - Professional and clear */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
        margin-bottom: 1.5rem;
        border-bottom: 2px solid var(--deep-red); /* Stronger separator for main tabs */
    }
    .stTabs [data-baseweb="tab"] {
        background-color: var(--soft-red); /* Light red for inactive tabs */
        border-radius: 8px 8px 0 0;
        padding: 10px 20px;
        color: var(--text-darker-red); /* Dark text on light red */
        font-weight: 500;
        transition: background-color 0.3s, color 0.3s;
        border: 1px solid var(--soft-red); /* Subtle border */
        border-bottom: none; /* No bottom border for inactive */
    }
    .stTabs [data-baseweb="tab"]:hover {
        background-color: var(--light-accent-red); /* Lighter red on hover */
        color: var(--text-color-dark); /* White text on hover */
    }
    .stTabs [aria-selected="true"] {
        background-color: var(--primary-red); /* Primary red for active tab */
        color: var(--text-color-dark); /* White text on active tab */
        border: 1px solid var(--primary-red); /* Border matches active color */
        border-bottom: none; /* No bottom border for active */
        font-weight: 600; /* Bolder for active tab */
    }

    /* Buttons - Clear action, professional look */
    .stButton>button {
        background-color: var(--primary-red);
        color: var(--text-color-dark); /* White text */
        border-radius: 8px;
        border: none;
        padding: 0.7em 1.2em;
        font-weight: 600;
        transition: background-color 0.3s, transform 0.2s;
        box-shadow: 0 2px 5px var(--shadow-color);
        cursor: pointer;
    }
    .stButton>button:hover {
        background-color: var(--dark-red); /* Darker red on hover */
        color: #FFF;
        transform: translateY(-2px);
    }
    .stButton>button:active {
        background-color: var(--dark-red);
        color: #FFF;
        transform: translateY(0);
        box-shadow: none;
    }

    /* File Uploader - Clear visual cue */
    .stFileUploader {
        border: 2px dashed var(--medium-gray); /* Dashed border for upload area */
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
        background-color: var(--card-bg);
        margin-top: 1rem;
        margin-bottom: 1rem;
    }

    /* Expander/Collapsible Panels - Clean and distinguishable */
    
    .expander_div {
        background-color: var(--soft-gray); /* Soft gray header */
        border-radius: 8px;
        padding: 10px 15px;
        font-weight: 500;
        color: var(--text-color-light);
        transition: background-color 0.3s;
        border: 1px solid var(--soft-gray);
        font-size: 1.2rem;
    }
    .expander_div:hover {
        background-color: #d0d0d0; /* Slightly darker gray on hover */
    }
    .expander_subdiv {
        background-color: var(--card-bg); /* White content area */
        border: 1px solid var(--soft-gray);
        border-top: none; /* No top border to blend with header */
        border-radius: 0 0 8px 8px;
        padding: 15px;
    }

    /* Metrics - Highlighted key information */
    [data-testid="stMetric"] {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1rem;
        box-shadow: 0 2px 5px var(--shadow-color);
        border: 1px solid var(--border-color);
        margin-bottom: 1rem;
    }
    [data-testid="stMetricLabel"] {
        color: var(--primary-red); /* Red label for importance */
        font-weight: 600;
        font-size: 0.9em;
    }
    [data-testid="stMetricValue"] {
        color: var(--text-darker-red); /* Darker red for metric values */
        font-size: 1.8em;
        font-weight: 700;
    }
    [data-testid="stMetricDelta"] {
        color: #4CAF50; /* Green for positive delta, adjust for negative if needed */
    }

    /* Text Input & Text Area - Clean borders, clear focus */
    .stTextInput>div>div>input, .stTextArea>div>div>textarea {
        border: 1px solid var(--medium-gray); /* Medium gray border */
        border-radius: 8px;
        padding: 8px 12px;
        transition: border-color 0.3s, box-shadow 0.3s;
        color: var(--text-color-light);
        background-color: var(--card-bg); /* Explicitly set background to white */
    }
    .stTextInput>div>div>input:focus, .stTextArea>div>div>textarea:focus {
        border-color: var(--primary-red); /* Red border on focus */
        box-shadow: 0 0 0 1px var(--primary-red); /* Red shadow on focus */
        outline: none; /* Remove default outline */
    }

    /* Info/Warning/Error Messages - Clear visual cues */
    .stAlert {
        border-radius: 8px;
        margin-bottom: 1rem;
    }
    .stAlert > div > div {
        /* Default background and text color for info */
        background-color: #e6f7ff; 
        color: #004085;
        border-left: 5px solid #007bff;
    }
    .stAlert.st-warning > div > div { /* Specific for warning */
        background-color: #fff3cd; 
        color: #664d03;
        border-left: 5px solid #ffc107;
    }
    .stAlert.st-success > div > div { /* Specific for success */
        background-color: #d4edda; 
        color: #0f5132;
        border-left: 5px solid #28a745;
    }
    .stAlert.st-error > div > div { /* Specific for error */
        background-color: #f8d7da; 
        color: #842029;
        border-left: 5px solid #dc3545;
    }

    

    /* General containers for visual separation - "Panels" */
    .stContainer {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1.5rem;
        margin-bottom: 1.5rem;
        box-shadow: 0 2px 5px var(--shadow-color);
        border: 1px solid var(--border-color);
    }

    /* Dataframes - Clean and readable tables */
    .stDataFrame {
        border: 1px solid var(--medium-gray);
        border-radius: 8px;
        overflow: hidden; /* Ensures border-radius applies to content */
    }
    .stDataFrame > div > div {
        background-color: var(--card-bg);
    }
    .stDataFrame th {
        background-color: var(--soft-gray);
        color: var(--text-color-light);
        font-weight: 600;
        border-bottom: 1px solid var(--medium-gray);
    }
    .stDataFrame td {
        color: var(--text-color-light);
        border-bottom: 1px solid var(--soft-gray);
    }
    .stDataFrame tr:nth-child(even) {
        background-color: #f5f5f5; /* Subtle stripe for readability */
    }

    
     /* sidebar style*/
    [data-testid="stSidebar"] {
        background-color: #0e1117;
        border-right: 1px solid var(--medium-gray);
        padding-top: 1rem;
    }
    
    /* Sidebar header*/
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h1,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h2,
    [data-testid="stSidebar"] [data-testid="stMarkdown"] h3 {
        color: var(--deep-red);
        font-weight: 600;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid var(--light-red);
        margin-bottom: 1.5rem;
    }
    
    /* Radio buttons on the sidebar container */
    [data-testid="stSidebar"] .stRadio {
        margin-bottom: 2rem;
    }
    
    /* Style for options buttons */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label > div {
        background-color: var(--sidebar-item-bg);
        border-radius: 8px;
        padding: 12px 15px;
        margin-bottom: 0.75rem;
        transition: all 0.2s ease;
        color: #FFF;
        font-weight: 500;
        border-left: 4px solid transparent;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    
    /* Hover state for options for radio button */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label:hover > div {
        background-color: var(--sidebar-item-hover);
        border-left: 4px solid var(--light-red);
        transform: translateX(3px);
    }
    
    /* Active state (selected) for radio buttom options */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] input[type="radio"]:checked + div {
        background-color: var(--sidebar-item-active);
        color: var(--text-color-dark);
        border-left: 4px solid var(--deep-red);
        box-shadow: 0 2px 8px rgba(229, 41, 32, 0.3);
        transform: translateX(5px);
    }
    
    /* Sidebar sections */
    [data-testid="stSidebar"] .sidebar-section {
        background-color: var(--card-bg);
        border-radius: 8px;
        padding: 1rem;
        margin-bottom: 1.5rem;
        border: 1px solid var(--medium-gray);
        box-shadow: 0 2px 5px rgba(0,0,0,0.05);
    }
    
    /* Section title sidebar */
    [data-testid="stSidebar"] .sidebar-section-title {
        color: var(--deep-red);
        font-weight: 600;
        font-size: 1rem;
        margin-bottom: 0.75rem;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid var(--soft-red);
    }
    
    /* Sidebar buttons */
    [data-testid="stSidebar"] .stButton > button {
        background-color: var(--soft-red);
        color: #FFF;
        border: 1px solid var(--light-red);
        border-radius: 6px;
        font-weight: 500;
        transition: all 0.2s ease;
        width: 100%;
        margin-bottom: 0.5rem;
    }
    
    /* Hover state for botones on sidebar */
    [data-testid="stSidebar"] .stButton > button:hover {
        background-color: var(--light-red);
        color: #FFF;
        border-color: var(--primary-red);
        transform: translateY(-1px);
        box-shadow: 0 2px 5px rgba(229, 41, 32, 0.2);
    }
    
    
    
    
    /* Sidebar Metrix*/
    [data-testid="stSidebar"] [data-testid="stMetric"] {
        background-color: var(--soft-red);
        padding: 0.75rem;
        border-radius: 6px;
        margin-bottom: 0.75rem;
    }
    
    [data-testid="stSidebar"] [data-testid="stMetric"] label {
        color: var(--deep-red) !important;
        font-weight: 600 !important;
    }
    
    [data-testid="stSidebar"] [data-testid="stMetric"] [data-testid="metric-container"] {
        color: var(--text-color-light) !important;
    }
    
    /* Dataframes - Clean and readable tables */
    .stDataFrame {
        border: 1px solid var(--medium-gray);
        border-radius: 8px;
        overflow: hidden; /* Ensures border-radius applies to content */
    }
    .stDataFrame > div > div {
        background-color: var(--card-bg);
    }
    .stDataFrame th {
        background-color: var(--soft-gray);
        color: var(--text-color-light);
        font-weight: 600;
        border-bottom: 1px solid var(--medium-gray);
    }
    .stDataFrame td {
        color: var(--text-color-light);
        border-bottom: 1px solid var(--soft-gray);
    }
    .stDataFrame tr:nth-child(even) {
        background-color: #f5f5f5; /* Subtle stripe for readability */
    }
    /* Main header style */
    .main-header {
        background: linear-gradient(135deg, var(--header-gradient-start), var(--header-gradient-end));
        border-radius: 12px;
        box-shadow: 0 8px 24px var(--header-shadow);
        margin-bottom: 2rem;
        padding: 2.5rem 2rem;
        border: 1px solid var(--header-border);
        position: relative;
        overflow: hidden;
    }
    
    /* Header decorations */
    .header-decoration {
        position: absolute;
        top: 0;
        right: 0;
        width: 100%;
        height: 100%;
        background-image: 
            linear-gradient(45deg, var(--text-darker-red) 25%, transparent 25%),
            linear-gradient(-45deg, var(--text-darker-red) 25%, transparent 25%);
        background-size: 12px 12px;
        background-position: 0 0, 6px 0;
        opacity: 0.03;
        z-index: 0;
        clip-path: polygon(70% 0, 100% 0, 100% 100%, 90% 100%);
    }

    .header-decoration-2 {
        position: absolute;
        bottom: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-image: 
            linear-gradient(45deg, var(--soft-red) 25%, transparent 25%),
            linear-gradient(-45deg, var(--soft-red) 25%, transparent 25%);
        background-size: 10px 10px;
        background-position: 0 0, 5px 0;
        opacity: 0.02;
        z-index: 0;
        clip-path: polygon(0 30%, 0 100%, 30% 100%);
    }
    
    /* Header content */
    .header-content {
        position: relative;
        z-index: 1;
        text-align: center;
    }
    
    /* MAin title */
    .app-title {
        font-size: 3rem;
        font-weight: 700;
        color: var(--deep-red);
        margin-bottom: 0.75rem;
        letter-spacing: -0.5px;
        text-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    
    /* subtitle */
    .app-description {
        font-size: 1.15rem;
        color: var(--text-color-light);
        max-width: 800px;
        margin: 0 auto;
        line-height: 1.6;
    }
</style>
""", unsafe_allow_html=True)


# ---------------- MAIN APPLICATION LAYOUT ----------------

#Header
html_content = """
<style>
        /* Primary Color Palette based on #e52920 */
        :root {
            --bg-light-red: #fdeeee;           /* Very light red/pink for backgrounds */
            --soft-red: #fbdad9;               /* Light red/pink for inactive elements, light card borders */
            --medium-light-red: #f7bcbb;       /* Slightly darker pink */
            --light-accent-red: #f28f8c;       /* Medium pink/light red for hover states */
            --primary-red: #ea534e;            /* Adjusted to #ea534e as per new palette, was #e52920 */
            --dark-accent-red: #e62e25;        /* Slightly darker primary red */
            --dark-red: #e41513;               /* Darker red for button hover/active */
            --very-dark-red: #be1818;          /* Even darker red */
            --deep-red: #9d1915;               /* Dark red/maroon for strong accents/borders */
            --text-darker-red: #53110e;        /* Very dark red for text/headers */
            
            --soft-gray: #e0e0e0;              /* Soft gray for general backgrounds/borders (kept as is) */
            --medium-gray: #cccccc;            /* Medium gray for borders/dividers (kept as is) */
            --dark-gray: #1a1a1a;              /* General dark text/header color (kept as is) */
            --text-color-light: #333333;       /* Dark text on light backgrounds (kept as is) */
            --text-color-dark: #ffffff;        /* White text on dark primary colors (kept as is) */
            --card-bg: #ffffff;                /* White background for cards/panels (kept as is) */
            --border-color: var(--soft-gray);  /* Default border color (kept as is) */
            --shadow-color: rgba(0, 0, 0, 0.1); /* Subtle shadow (kept as is) */
        }

        /* Headers - Ensure good contrast */
        h1, h2, h3, h4, h5, h6 {
            color: var(--text-darker-red); /* Using the darkest red for headers */
            font-weight: 600;
        }
        h1 { 
            font-size: 2.5em; 
            color: var(--primary-red); /* Keep H1 primary red */
            margin-bottom: 0.5em; 
        } 
        h2 { font-size: 2em; margin-bottom: 0.5em; }
        h3 { font-size: 1.5em; margin-bottom: 0.5em; }


        /* Existing styles adapted to new variables */
        body {
            font-family: 'Poppins', sans-serif;
            margin: 0;
        }
        
        .main-header {
            padding: 2rem 1rem;
            background-color: var(--primary-red);
            border-radius: 2rem;
            text-align: center;
        }
        .app-title {
            color: var(--medium-light-red);
            margin-bottom: 0.5rem;
        }
        .app-description, .secondary-line {
            font-size: 1.4rem;
            color: #000000;
            margin: 0.3rem 0;
        }
        .rotating-quote {
            font-size: 1.3rem;
            font-style: italic;
            margin-top: 1.2rem;
            color: var(--text-darker-red);
            min-height: 2rem;
            transition: opacity 0.3s ease-out, transform 0.3s ease-out;
            opacity: 1; 
            transform: translateY(0); 
        }

        /* --- RESPONSIVE STYLES --- */
        @media (max-width: 768px) {
            h1 {
                font-size: 2em;
            }
            .main-header {
                padding: 1.5rem 0.5rem;
            }
            .app-description, .secondary-line {
                font-size: 1.2rem;
            }
            .rotating-quote {
                font-size: 1rem;
                min-height: 1.5rem;
            }
        }
        
        @media (max-width: 480px) {
            h1 {
                font-size: 1.5em;
            }
            .main-header {
                border-radius: 1rem;
            }
            .app-description, .secondary-line {
                font-size: 1rem;
            }
        }
    </style>
</head>
<body>
    <div class="main-header">
        <div class="header-content">
            <h1 class="app-title">TADPOLE: Build your own RNA switch.</h1>
            <p class="app-description">Get ON/OFF systems that reprogram translation through shape, not sequence.</p>
            <p class="secondary-line">Input your sequence. TADPOLE finds the path to control.</p>
            <p class="rotating-quote" id="quote"></p>
        </div>
    </div>

    <script>
        const quotes = [
            "Forget promoters. Think structure.<br>The future of regulation folds differently.",
            "What if translation was the new frontier of control?<br>Welcome to the RNA switch revolution.",
            "Biology has always known how to self-regulate.<br>We’re just learning to speak its language.",
            "Fast, reversible, translation-level regulation?<br>We’re already skipping limits. Join us, beyond the finish line."
        ];
        let quoteIndex = 0;
        const quoteElement = document.getElementById("quote");

        function rotateQuote() {
            quoteElement.style.opacity = 0;
            quoteElement.style.transform = 'translateY(10px)'; 

            setTimeout(() => {
                quoteElement.innerHTML = quotes[quoteIndex];
                quoteIndex = (quoteIndex + 1) % quotes.length;
                
                quoteElement.style.transform = 'translateY(-10px)'; 

                setTimeout(() => {
                    quoteElement.style.opacity = 1;
                    quoteElement.style.transform = 'translateY(0)'; 
                }, 10); 

            }, 300); 
        }

        document.addEventListener('DOMContentLoaded', () => {
            rotateQuote();
            setInterval(rotateQuote, 6000); 
        });
    </script>

"""

st.components.v1.html(html_content, height=350)

# Sidebar for navegation
with st.sidebar:
    st.image("tadpole_package/images/logo.png")
    st.header("TADPOLE")
    selected_main_tab = st.radio(
        "Navegation",
        ["Tools", "Help"],
        index=1 # default: "Model System"
    )

# Main content (acoding to sidebar selection)
if selected_main_tab == "Help":
    app_help.help_tab()
elif selected_main_tab == "Tools":

    tools_sub_tabs = st.tabs(["Model System", "Structural RNA element"])

    with tools_sub_tabs[0]: # "Model system"
        app_switch_design.switch_design() 
    with tools_sub_tabs[1]: # "structural RNA element"
        app_SRE.structural_rna_element_tab() 

   


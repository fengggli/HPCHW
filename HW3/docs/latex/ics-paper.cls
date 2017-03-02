\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ics-paper}[2016/12/28 ICS paper template]

% Load the underlying article class
\LoadClass[letterpaper,10pt,twocolumn]{article}

% Page margins
\RequirePackage{geometry}
\geometry{lmargin=0.75in,rmargin=0.75in,tmargin=1in,bmargin=1in}
\setlength{\columnsep}{0.25in}

% Font
\usepackage{times}

% Title
\RequirePackage{titling}
\setlength{\droptitle}{-3.5em}

% Abstract
\RequirePackage{abstract}
\renewcommand{\abstractnamefont}{\large\bfseries\sffamily}
\renewcommand{\abstracttextfont}{\normalfont}
\setlength{\absleftindent}{0in}
\setlength{\absrightindent}{0in}

% Macro
\RequirePackage{xspace}

% Section styles
\RequirePackage{titlesec}
\titleformat{\section}{\large\bfseries\sffamily}{\thesection. }{0pt}{}
\titleformat{\subsection}{\bfseries\sffamily}{\thesubsection. }{0pt}{}
\titleformat{\subsubsection}{\sffamily}{\thesubsubsection. }{0pt}{}

% Paragraphs
\setlength{\parindent}{1em}

% Captions
\RequirePackage[font={sf,small,bf}]{caption}
\RequirePackage{subcaption}
\setlength{\abovecaptionskip}{0.5em}

% Tables
\RequirePackage{array}
\RequirePackage{multirow}
\RequirePackage{booktabs}

% Figures
\RequirePackage{graphicx}
\RequirePackage{dblfloatfix}

% Footnotes
\RequirePackage{footmisc}
\setlength{\footnotemargin}{2pt}

% Bibliography
\RequirePackage[numbers,sort]{natbib}

% URL package
\RequirePackage{url}
\usepackage{hyperref}

\endinput

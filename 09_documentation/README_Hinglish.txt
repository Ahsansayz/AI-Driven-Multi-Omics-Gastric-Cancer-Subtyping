================================================================================
       GASTRIC CANCER MUTATIONAL SIGNATURE ANALYSIS — HINGLISH GUIDE
       (Ye poora project kya hai, kyun hai, aur kaise kiya — sab samjho)
================================================================================


SABSE PEHLE — YE PROJECT HAI KYA?
══════════════════════════════════

Dekho bhai, ye project Gastric Cancer (stomach ka cancer) ke baare mein hai.

Cancer hota kaise hai? Tumhare DNA mein mutations aate hain — matlab kuch 
letters (A, T, G, C) galat ho jaate hain. Ab har cancer mein alag tarah ke 
mutations hote hain. Kuch mutations aging ki wajah se hote hain, kuch 
smoking se, kuch body ke andar ke repair system kharab hone se.

In mutation patterns ko "MUTATIONAL SIGNATURES" kehte hain.

TCGA (The Cancer Genome Atlas) ne bataya ki gastric cancer ke 4-5 subtypes 
hote hain:
  - CIN  → Chromosomally Instable (sabse common, ~58%)
  - MSI  → Microsatellite Instable (DNA repair kharab, ~20%)
  - GS   → Genomically Stable (stable genome, diffuse type, ~12%)
  - EBV  → Epstein-Barr Virus positive (virus wala, ~8%)
  - POLE → Ultra-mutated (bahut zyada mutations, rare ~2%)

Humne kya kiya? Humne AI/ML se sirf mutation data dekh ke predict kiya ki 
konsa patient kis subtype mein aata hai. Ye important hai kyunki alag 
subtype ka alag treatment hota hai (jaise MSI patients ko immunotherapy 
acchi lagti hai).


STEP 0: SETUP — PEHLE KYA TAYYARI KI
══════════════════════════════════════

(a) DATA KAHAN SE LIYA?
    - Website: https://portal.gdc.cancer.gov
    - Project: TCGA-STAD (Stomach Adenocarcinoma)
    - Kya download kiya: 434 MAF files
    - MAF file kya hai? Ye ek badi table hai jismein ek patient ke saare 
      somatic mutations listed hain — konsa gene, konsa chromosome, kya 
      change hua, sab kuch.
    - Har file ek patient ki hai.

(b) FOLDER STRUCTURE BANAYA:
    mkdir -p /home/param/Downloads/Project/maf_files
    mkdir -p /home/param/Downloads/Project/output
    mkdir -p /home/param/Downloads/Project/data
    
    Saari MAF files "maf_files/" folder mein daali.

(c) PYTHON PACKAGES INSTALL KIYE:
    pip install pandas numpy scikit-learn xgboost matplotlib seaborn \
      tqdm shap imbalanced-learn lifelines umap-learn
    pip install SigProfilerAssignment
    pip install python-docx

    Ye sab tools hain — pandas data ke liye, scikit-learn ML ke liye,
    shap model samajhne ke liye, lifelines survival analysis ke liye.


STEP 1: SBS96 MATRIX BANANA — MUTATION KA FINGERPRINT
══════════════════════════════════════════════════════

SAMJHO KYA HUA:
    Dekho, agar tumhare DNA mein "C" letter "T" mein badal gaya, to ye ek 
    C>T mutation hai. Lekin SIRF ye jaanna kaafi nahi — humein ye bhi 
    dekhna hai ki us C ke dono taraf kya tha.

    Jaise: A_C_G mein C>T hua → ye hai "A[C>T]G"
           T_C_A mein C>T hua → ye hai "T[C>T]A"

    Ye dono ALAG mutation types hain! Kyunki surrounding bases matter 
    karte hain — wo batate hain ki KONSI biological process ne ye 
    mutation cause kiya.

    Total kitne combinations ban sakte hain?
    - 6 substitution types: C>A, C>G, C>T, T>A, T>C, T>G
    - Har ek ke 16 trinucleotide contexts (4 left × 4 right bases)
    - Total: 6 × 16 = 96 channels

    Isko kehte hain "96-channel SBS trinucleotide context matrix"
    — ye GOLD STANDARD hai mutational signature analysis ka.

KYA KIYA:
    - Har MAF file se SNPs (single base changes) nikale
    - CONTEXT column se trinucleotide context extract kiya
    - Pyrimidine convention apply kiya (agar ref base A ya G hai to 
      complement le lo — C ya T banao)
    - 96 channels mein count kiya

COMMAND:
    python step1_build_sbs96_matrix.py

RESULT:
    - 431 samples × 96 mutation types ki matrix bani
    - Total 1,54,534 mutations capture hue
    - Sabse common mutation: G[C>T]G (20,420 baar) — ye aging ka sign hai
    - Save hua: output/sbs96_matrix.csv


STEP 2: SIGNATURES NIKALNA — NMF SE
════════════════════════════════════

SAMJHO KYA HUA:
    Ab humne socha — in 431 patients mein konse "patterns" hain? Jaise 
    agar 50 patients mein same tarah ka mutation profile hai, to shayad 
    unme same mutagenic process chal raha hai.

    Isko nikalne ke liye NMF (Non-negative Matrix Factorization) use kiya.
    
    NMF kya karta hai? Simple bhasha mein:
    Tumhari badi matrix (431 patients × 96 channels) ko 2 chhoti matrices 
    mein tod deta hai:
      - W matrix: Har patient mein kitna konsa signature active hai
      - H matrix: Har signature ka 96-channel profile kya dikhta hai

    Lekin sawaal ye hai: KITNE signatures hain? 2? 5? 10? 15?

    Iske liye humne "rank selection" kiya — NMF ko 2 se 12 tak test kiya,
    har rank pe 50 baar run kiya, aur dekha:
    1. Reconstruction error kitna kam ho raha hai (kam = better)
    2. Stability kitni hai (agar 50 baar run karne pe same answer aaye = stable)

    Result: k=9 optimal nikla (9 signatures).

COMMAND:
    python step2_extract_signatures.py

RESULT:
    - 9 de novo signatures extract hue
    - Saari ranks pe stability 1.0 thi (excellent!)
    - Reconstruction error: 0.079 (matlab 92% data explain ho raha hai)
    - Plots save hue: rank_selection.png, signature_profiles_96ch.png


STEP 3: COSMIC SIGNATURES SE MATCH KARNA
═════════════════════════════════════════

SAMJHO KYA HUA:
    Ab humne apne 9 de novo signatures nikale — lekin inke koi naam nahi 
    hain. Scientists ne already 86 known signatures catalogue kiye hain 
    (COSMIC database) — jaise SBS1 = aging, SBS6 = DNA repair defect, etc.

    To humne kya kiya? Har patient ka 96-channel profile liya aur COSMIC 
    ke 86 reference signatures se match kiya — "ye patient mein SBS1 
    kitna hai, SBS6 kitna hai, SBS15 kitna hai..." ye sab nikala.

    Method: NNLS (Non-Negative Least Squares)
    Simple mein: Ye ek math formula hai jo best possible combination 
    dhundta hai COSMIC signatures ka, taaki tumhare patient ka profile 
    ban sake. Condition: sab weights >= 0 hone chahiye (negative nahi).

    Quality check: Cosine Similarity — original profile aur reconstruction 
    kitne similar hain? 1.0 = perfect match.

PEHLE COSMIC REFERENCE FILE LENI PADI:
    pip install SigProfilerAssignment
    # Iske andar COSMIC v3.4 reference bundled hai:
    # 96 channels × 86 signatures

COMMAND:
    python step3_cosmic_assignment.py

RESULT:
    - 431 samples fitted to 86 COSMIC signatures
    - 83 signatures active mile (3 signatures kisi mein bhi nahi the)
    - Mean cosine similarity: 0.933 (bahut accha! 0.90+ excellent maana jata hai)
    - 348 samples (80.7%) ne 0.90+ achieve kiya
    - Top signatures: SBS15 (DNA repair defect), SBS1 (aging), SBS6 (repair)
    - Save hua: output/cosmic_assignment/cosmic_activities.csv


STEP 4: CLINICAL DATA DOWNLOAD KARNA
═════════════════════════════════════

SAMJHO KYA HUA:
    Ab tak humne sirf mutations dekhe. Lekin ML model train karne ke liye 
    humein LABELS chahiye — matlab ye jaanna ki "ye patient CIN hai, ye 
    MSI hai, ye EBV hai."

    Ye labels TCGA ne apni 2014 Nature paper mein publish kiye the.
    Humne 2 jagah se data liya:

    (a) GDC API (Government cancer database):
        - Gender, age, vital status, days to death
        - POST request bheja (kyunki GET mein URL bahut lamba ho jata tha)
        - 443 patients ka data mila

    (b) cBioPortal API (Open cancer genomics portal):
        - Molecular subtype labels (CIN/MSI/GS/EBV/POLE)
        - 383 patients ke official labels mile
        
    (c) TMB (Tumor Mutational Burden) calculate kiya:
        TMB = mutations / 30 Mb (exome size)
        - Ye batata hai ki tumor kitna "mutated" hai
        - MSI patients ka TMB bahut high hota hai

    Sab merge karke 375 patients ke paas label + mutation data dono mil gaya.

COMMAND:
    python step4_get_clinical_data.py

RESULT:
    Subtype distribution:
      CIN:  219 patients (58.4%) — sabse common
      MSI:   73 patients (19.5%) — high mutation rate
      GS:    46 patients (12.3%) — stable genome
      EBV:   30 patients (8.0%)  — virus positive
      POLE:   7 patients (1.9%)  — ultra-rare
    
    Gender: 276 male, 155 female
    TMB: mean=10.6, median=3.2 mutations/Mb


STEP 5: FEATURE ENGINEERING — ML KE LIYE DATA TAYYAR KARNA
═══════════════════════════════════════════════════════════

SAMJHO KYA HUA:
    ML model ko sirf raw COSMIC activities dena kaafi nahi — humne kuch 
    EXTRA SMART FEATURES bhi banaye jo biology se inspired hain:

    179 total features bane:
    ├── 83 raw COSMIC signature activities (absolute counts)
    ├── 83 normalized proportions (percentage — 0 to 1)
    └── 13 engineered features:
        ├── MSI_sig_burden       → SBS6+SBS14+SBS15+SBS20+SBS21+SBS26+SBS44
        ├── APOBEC_burden        → SBS2+SBS13 (APOBEC enzyme activity)
        ├── clock_burden         → SBS1+SBS5 (aging signatures)
        ├── HRD_burden           → SBS3 (DNA repair by homologous recomb.)
        ├── SBS17_burden         → SBS17a+SBS17b (gastric-specific)
        ├── SBS1/SBS5 ratio      → aging vs background ratio
        ├── TMB                  → mutations per megabase
        ├── log_TMB              → log scale TMB
        ├── log_total_activity   → total mutation activity (log)
        ├── dominant_sig_prop    → sabse bada signature kitna dominant hai
        └── n_active_sigs        → kitne signatures active hain (>5%)

    YE KYUN IMPORTANT HAI?
    Raw data se model ko patterns dhundne mein mushkil hoti hai. Jab hum 
    biology ka knowledge use karke features banate hain (jaise "saare MMR 
    signatures ka sum"), to model ko bahut easy ho jata hai.

COMMAND:
    python step5_feature_engineering.py

RESULT:
    - 375 labeled samples × 179 features ki matrix bani
    - Save hua: output/ml_features.csv, output/ml_labels.csv


STEP 6: ML CLASSIFICATION — AI MODEL TRAINING
══════════════════════════════════════════════

SAMJHO KYA HUA:
    Ab aaya main step! 5 alag ML models train kiye:

    1. Random Forest — bahut saare decision trees ka ensemble
    2. XGBoost — gradient boosted trees (competitions mein champion)
    3. SVM — Support Vector Machine (boundary dhundta hai classes ke beech)
    4. MLP — Neural Network (deep learning lite)
    5. Gradient Boosting — sequential tree building

    Har model ke liye:
    - StandardScaler se data normalize kiya (sab features same scale pe)
    - GridSearchCV se best hyperparameters dhundhe (inner 3-fold CV)
    - 5-fold Stratified Cross-Validation se evaluate kiya
      (data ko 5 parts mein baanto, 4 pe train, 1 pe test, 5 baar rotate)

    METRICS KYA DEKHE:
    - Accuracy: Kitne % sahi predict kiye
    - F1 Score: Precision aur Recall ka balance (imbalanced data ke liye best)
    - MCC: Matthews Correlation Coefficient (-1 se +1, best overall metric)
    - AUC-ROC: ROC curve ke neeche ka area (1.0 = perfect)

COMMAND:
    python step6_ml_classification.py

RESULTS:
    ┌───────────────────┬──────────┬────────┬───────┬─────────┐
    │ Model             │ Accuracy │ F1     │ MCC   │ AUC-ROC │
    ├───────────────────┼──────────┼────────┼───────┼─────────┤
    │ XGBoost 🏆 BEST   │ 81.1%    │ 0.765  │ 0.670 │ 0.870   │
    │ Random Forest     │ 80.3%    │ 0.744  │ 0.660 │ 0.873   │
    │ Gradient Boosting │ 80.0%    │ 0.761  │ 0.648 │ 0.844   │
    │ MLP Neural Net    │ 74.4%    │ 0.666  │ 0.539 │ 0.717   │
    │ SVM (RBF)         │ 72.0%    │ 0.673  │ 0.492 │ 0.812   │
    └───────────────────┴──────────┴────────┴───────┴─────────┘

    Per-subtype (XGBoost):
    - MSI:  100% correct!! (73/73) — kyunki MSI ka mutation pattern BAHUT unique hai
    - CIN:  97.3% correct (213/219) — majority class, accha perform kiya
    - GS:   32.6% correct (15/46) — low mutations, hard to distinguish
    - POLE: 28.6% correct (2/7) — sirf 7 samples the, too few
    - EBV:  3.3% correct (1/30) — small sample + overlapping features

    TOP 5 IMPORTANT FEATURES (model ne sabse zyada kin features pe rely kiya):
    1. SBS54  — MMR related signature
    2. SBS15  — DNA mismatch repair deficiency
    3. TMB    — Total mutation burden
    4. MSI_sig_burden_prop — Humara engineered feature!
    5. SBS23  — Unknown etiology


STEP 7: VISUALIZATION — PLOTS BANANA
═════════════════════════════════════

SAMJHO KYA HUA:
    Research paper/dissertation ke liye publication-quality figures chahiye.
    Humne 6 types ke plots banaye:

    1. SHAP PLOTS — "Model ne YE decision KYUN liya?"
       SHAP values batate hain ki har feature ne prediction ko kitna 
       aur kis direction mein push kiya. Jaise: "High SBS15 value ne 
       is patient ko MSI ki taraf push kiya."

    2. t-SNE + UMAP — "Kya subtypes alag clusters banate hain?"
       179-dimensional data ko 2D mein project kiya. Agar same subtype 
       ke patients ek jagah cluster karte hain = good separation.
       Result: MSI ka cluster CLEARLY alag dikh raha tha!

    3. SIGNATURE-SUBTYPE HEATMAP — "Konse signature konse subtype mein zyada?"
       Matrix dikhata hai ki MSI mein SBS15/SBS6 high hai, CIN mein 
       SBS17b high hai, etc.

    4. KAPLAN-MEIER SURVIVAL CURVES — "Konsa subtype zyada survive karta hai?"
       Overall Survival curves by subtype — kaunse patients zyada 
       time tak jeete hain.

    5. TMB DISTRIBUTION — "Konse subtype mein zyada mutations?"
       Box plots dikhate hain ki MSI/POLE mein TMB bahut high hai.

COMMAND:
    python step7_visualize.py

RESULT:
    6 PNG figures save hue output/figures/ mein


STEP 8: HTML REPORT GENERATE KARNA
═══════════════════════════════════

SAMJHO KYA HUA:
    Saare results, plots, tables ek professional dark-themed HTML report 
    mein compile kiye. Saari images base64 mein embed hain (matlab report 
    ek single file hai, koi external dependency nahi).

COMMAND:
    python step8_generate_report.py

RESULT:
    - output/report.html (5.4 MB)
    - Browser mein kholne ka command: xdg-open output/report.html


BONUS: RESEARCH PAPER (DOCX) GENERATE KARNA
════════════════════════════════════════════

    Ek full research paper DOCX format mein bhi generate kiya:
    - Title page, Abstract, Introduction, Methods, Results, Discussion
    - 16 figures embedded, 5 tables, 12 references
    - ~4,300 words, ~30 pages

COMMAND:
    python generate_paper.py

RESULT:
    output/Research_Paper_Gastric_Cancer_Mutational_Signatures.docx (3.4 MB)


POORA PIPELINE EK SAATH CHALANA HO TO:
═══════════════════════════════════════

    # Activate conda environment
    source /home/param/miniconda3/etc/profile.d/conda.sh
    conda activate base

    # Go to project folder
    cd /home/param/Downloads/Project

    # Run everything
    bash run_pipeline.sh

    # Ya sirf specific steps
    bash run_pipeline.sh 1       # sirf step 1
    bash run_pipeline.sh 6       # sirf ML training
    bash run_pipeline.sh 3-6     # step 3 se 6 tak


AGAR KISI KO SIMPLE MEIN SAMJHANA HO:
══════════════════════════════════════

    "Humne 431 stomach cancer patients ka mutation data liya (TCGA se).
     Usse 96 tarah ke mutation patterns (trinucleotide context) count kiye.
     Phir NMF se 9 mutational signatures nikale.
     In signatures ko known COSMIC signatures se match kiya.
     Phir official cancer subtypes (CIN/MSI/GS/EBV) ke labels download kiye.
     Signatures + TMB + engineered features mila ke 179 features banaye.
     5 ML models train kiye — XGBoost ne 81% accuracy se subtypes predict 
     kar diye, aur MSI subtype 100% sahi classify hua.
     SHAP analysis se samjha ki model ne ye decisions kyun liye.
     Survival analysis se confirm kiya ki subtypes clinically relevant hain.
     Sab kuch ek HTML report aur research paper mein compile kar diya."


================================================================================
                              — THE END —
================================================================================

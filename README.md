# 🧬 Genetic Risk Analyzer

A comprehensive C++ tool for analyzing genetic data to assess potential health risks based on known variants, age-specific penetrance, and polygenic risk scores.

## 📋 Overview

The Genetic Risk Analyzer is designed to help researchers and individuals understand their genetic predisposition to various diseases by analyzing genomic data from multiple file formats. It combines variant-specific risk assessments with polygenic risk score (PRS) calculations to provide a holistic view of genetic health risks.

## ✨ Features

- **Multiple Input Formats**: Supports VCF files, FASTA files, and manual DNA sequence input
- **Age-Specific Risk Assessment**: Calculates disease penetrance based on age-specific data
- **Polygenic Risk Scores**: Computes PRS for common complex diseases
- **Environmental & Risk Factors**: Incorporates family history, lifestyle factors, and environmental influences
- **Quality Control**: Validates genotype calls with quality scores and read depth metrics
- **Visual Reports**: Generates risk assessment plots with confidence intervals using Gnuplot
- **ACMG Classification Support**: Handles variant classifications (Pathogenic, Likely Pathogenic, VUS, etc.)

## 🛠️ Requirements

### Dependencies

- **C++ Compiler**: C++11 or later (g++, clang++)
- **nlohmann/json**: JSON parsing library ([GitHub](https://github.com/nlohmann/json))
- **Gnuplot**: For generating visualization plots

### System Requirements

- Linux, macOS, or Windows (with appropriate compiler)
- Minimum 2GB RAM
- Terminal/Command line access

## 📦 Installation

1. **Clone the repository**:
```bash
   git clone https://github.com/yourusername/genetic-risk-analyzer.git
   cd genetic-risk-analyzer
```

2. **Install nlohmann/json library**:
   
   Download the single-header file:
```bash
   wget https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp
```
   
   Or install via package manager:
```bash
   # Ubuntu/Debian
   sudo apt-get install nlohmann-json3-dev
   
   # macOS
   brew install nlohmann-json
```

3. **Install Gnuplot** (for visualization):
```bash
   # Ubuntu/Debian
   sudo apt-get install gnuplot
   
   # macOS
   brew install gnuplot
   
   # Windows
   # Download from http://www.gnuplot.info/
```

4. **Compile the program**:
```bash
   g++ -std=c++11 -o genetic_analyzer main.cpp -lpthread
```

## 🚀 Usage

### Basic Usage

Run the program:
```bash
./genetic_analyzer
```

The program will guide you through an interactive menu where you can:

1. **Choose input type**:
   - VCF file (Variant Call Format)
   - FASTA file
   - Manual DNA sequence

2. **Enter your age** for age-specific risk calculations

3. **Provide additional factors**:
   - Family history of diseases
   - Smoking status
   - Other environmental factors

### Input File Formats

#### VCF File Example
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    12345   .       A       G       30      PASS    .       GT:DP   0/1:25
chr2    67890   .       C       T       35      PASS    .       GT:DP   1/1:30
```

#### FASTA File Example
```
>Sample_Sequence
ATGCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
```

#### disease_data.json Structure
```json
[
  {
    "chromosome": "chr17",
    "position": 43044295,
    "ref": "G",
    "alt": "A",
    "disease": "Breast Cancer",
    "description": "BRCA1 pathogenic variant",
    "isAgeSpecific": true,
    "agePenetrance": [
      {"age": 30, "penetrance": 0.03},
      {"age": 40, "penetrance": 0.19},
      {"age": 50, "penetrance": 0.50}
    ],
    "alleleFrequency": 0.0001,
    "geneImpact": "High",
    "classification": "Pathogenic",
    "oddsRatio": 5.0
  }
]
```

## 📊 Output

The program generates:

1. **Console Report**: Detailed text-based risk assessment showing:
   - Detected variants and their genotypes
   - Disease-specific risk percentages
   - Polygenic risk scores for common diseases
   - Quality warnings for low-confidence calls

2. **Visual Plot** (`risk_plot.png`): 
   - Bar chart showing risk percentages with confidence intervals
   - Automatically opened after analysis (on supported systems)

3. **Data File** (`plot_data.txt`):
   - Tab-delimited file with risk data for further analysis

## ⚠️ Important Disclaimers

- **Not for Clinical Use**: This tool is for research and educational purposes only
- **Consult Healthcare Professionals**: Always discuss genetic test results with qualified healthcare providers
- **Data Privacy**: Handle genetic data responsibly and ensure proper security measures
- **Accuracy Limitations**: Results depend on the quality and completeness of input data and disease databases

## 🔬 Technical Details

### Risk Calculation Methodology

1. **Base Penetrance**: Retrieved from age-specific lookup tables
2. **Genotype Adjustment**: Modified based on heterozygous (0/1) or homozygous (1/1) status
3. **Factor Modulation**: Adjusted using environmental and risk factor multipliers
4. **PRS Calculation**: Sum of log-odds ratios across multiple variants

### Quality Control

- Minimum variant quality score: 20
- Minimum read depth: 10
- Warnings displayed for low-quality calls

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, fork the repository, and create pull requests.

### Areas for Contribution

- Additional disease databases
- Support for more file formats
- Improved visualization options
- Population-specific frequency databases
- Enhanced PRS algorithms

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **nlohmann/json**: Excellent JSON library for C++
- **Gnuplot**: Powerful plotting utility
- **ACMG Guidelines**: Framework for variant classification
- The genomics research community for open data standards

## 📧 Contact

For questions, suggestions, or collaborations, please open an issue on GitHub or contact [your-email@example.com].

## 🔗 Resources

- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [ACMG Variant Classification Guidelines](https://www.acmg.net/)
- [ClinVar Database](https://www.ncbi.nlm.nih.gov/clinvar/)
- [gnomAD Population Frequencies](https://gnomad.broadinstitute.org/)

---

**Remember**: Genetic information is complex and probabilistic. This tool provides risk estimates based on current scientific knowledge, which is constantly evolving. Always seek professional genetic counseling for personalized interpretation.

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <limits>
#include <iomanip>  // For better formatting
#include "json.hpp" // Include the JSON library
#include <fstream>  // For file handling
#include <algorithm> // For sorting and searching
#include <cmath>     // For mathematical functions
#include <regex>     // For regular expressions
#include <map>       // For ordered maps
#include <fcntl.h>    // For file control options
#include <thread>     // For multi-threading
#include <mutex>      // For thread synchronization
#include <future>     // For async processing

using json = nlohmann::json;
using namespace std;

// Expanded structure to hold variant data
struct Variant {
    string chromosome;
    int position;
    string ref;
    string alt;
    string disease;
    string description;
    bool isAgeSpecific;
    vector<pair<int, float>> agePenetrance;
    
    // New fields
    float alleleFrequency;
    string geneImpact;
    string classification; // ACMG guidelines: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
    map<string, float> populationFrequencies; // Population-specific allele frequencies
    map<string, float> environmentalFactors; // Environmental factors that modulate risk
    map<string, float> riskFactors; // Other risk factors (e.g., family history)
    vector<float> confidenceIntervals; // Confidence intervals for penetrance
    float oddsRatio; // For polygenic risk score calculation
};

// Structure to hold genotype data from VCF
struct Genotype {
    string chromosome;
    int position;
    string ref;
    string alt;
    string genotype; // 0/0 (homozygous reference), 0/1 (heterozygous), 1/1 (homozygous alternate)
    float quality;
    int readDepth;
};

// Function to parse JSON file
vector<Variant> parseJSON(const string &filename) {
    vector<Variant> variants;
    ifstream file(filename);
    json jsonData;

    if (!file.is_open()) {
        cerr << "Error: Could not open the JSON file: " << filename << endl;
        return variants;
    }

    try {
        file >> jsonData;  // Parse the JSON data
    } catch (json::parse_error &e) {
        cerr << "Error: Could not parse the JSON file. " << e.what() << endl;
        return variants;
    }

    for (const auto &entry : jsonData) {
        try {
            Variant variant = {
                entry.at("chromosome").get<string>(),
                entry.at("position").get<int>(),
                entry.at("ref").get<string>(),
                entry.at("alt").get<string>(),
                entry.at("disease").get<string>(),
                entry.at("description").get<string>(),
                entry.at("isAgeSpecific").get<bool>(),
                {}
            };

            // Extract age penetrance data
            for (const auto &ageEntry : entry.at("agePenetrance")) {
                int age = ageEntry.at("age").get<int>();
                float penetrance = ageEntry.at("penetrance").get<float>();
                variant.agePenetrance.emplace_back(age, penetrance);
            }

            // Extract new fields if they exist
            if (entry.contains("alleleFrequency")) {
                variant.alleleFrequency = entry.at("alleleFrequency").get<float>();
            }
            if (entry.contains("geneImpact")) {
                variant.geneImpact = entry.at("geneImpact").get<string>();
            }
            if (entry.contains("classification")) {
                variant.classification = entry.at("classification").get<string>();
            }
            if (entry.contains("populationFrequencies")) {
                for (const auto &popEntry : entry.at("populationFrequencies").items()) {
                    variant.populationFrequencies[popEntry.key()] = popEntry.value().get<float>();
                }
            }
            if (entry.contains("environmentalFactors")) {
                for (const auto &envEntry : entry.at("environmentalFactors").items()) {
                    variant.environmentalFactors[envEntry.key()] = envEntry.value().get<float>();
                }
            }
            if (entry.contains("riskFactors")) {
                for (const auto &riskEntry : entry.at("riskFactors").items()) {
                    variant.riskFactors[riskEntry.key()] = riskEntry.value().get<float>();
                }
            }
            if (entry.contains("confidenceIntervals")) {
                for (const auto &ciEntry : entry.at("confidenceIntervals")) {
                    variant.confidenceIntervals.push_back(ciEntry.get<float>());
                }
            }
            if (entry.contains("oddsRatio")) {
                variant.oddsRatio = entry.at("oddsRatio").get<float>();
            }

            variants.push_back(variant);
        } catch (json::exception &e) {
            cerr << "Error: Missing data in JSON entry. " << e.what() << endl;
            continue;  // Skip the invalid entry
        }
    }

    file.close();
    return variants;
}

// Function to parse VCF file
vector<Genotype> parseVCF(const string &filename) {
    vector<Genotype> genotypes;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error: Could not open the VCF file: " << filename << endl;
        return genotypes;
    }

    // Skip header lines
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse VCF line
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, sample;
        int pos;
        
        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format >> sample)) {
            cerr << "Error: Invalid VCF line: " << line << endl;
            continue;
        }

        // Parse genotype
        string genotype;
        size_t gtIndex = format.find("GT");
        if (gtIndex != string::npos) {
            // Extract genotype from the sample field
            vector<string> formatFields;
            istringstream formatStream(format);
            string field;
            while (getline(formatStream, field, ':')) {
                formatFields.push_back(field);
            }

            vector<string> sampleFields;
            istringstream sampleStream(sample);
            while (getline(sampleStream, field, ':')) {
                sampleFields.push_back(field);
            }

            if (gtIndex < sampleFields.size()) {
                genotype = sampleFields[gtIndex];
            }
        }

        // Parse quality and read depth
        float quality = stof(qual);
        int readDepth = 0;
        
        size_t dpIndex = format.find("DP");
        if (dpIndex != string::npos) {
            // Extract read depth from the sample field
            vector<string> formatFields;
            istringstream formatStream(format);
            string field;
            while (getline(formatStream, field, ':')) {
                formatFields.push_back(field);
            }

            vector<string> sampleFields;
            istringstream sampleStream(sample);
            while (getline(sampleStream, field, ':')) {
                sampleFields.push_back(field);
            }

            if (dpIndex < sampleFields.size()) {
                readDepth = stoi(sampleFields[dpIndex]);
            }
        }

        // Add genotype to vector
        genotypes.push_back({chrom, pos, ref, alt, genotype, quality, readDepth});
    }

    file.close();
    return genotypes;
}

// Function to parse FASTA file
string parseFASTA(const string &filename) {
    ifstream file(filename);
    string line;
    string sequence;

    if (!file.is_open()) {
        cerr << "Error: Could not open the FASTA file: " << filename << endl;
        return sequence;
    }

    // Skip header line
    if (getline(file, line) && line[0] == '>') {
        // Read sequence lines
        while (getline(file, line)) {
            sequence += line;
        }
    }

    file.close();
    return sequence;
}

// Function to calculate penetrance based on age, genotype, and other factors
float calculatePenetrance(const Variant &variant, int age, const string &genotype, 
                          const map<string, string> &userFactors) {
    float basePenetrance = 0.0;

    // Get base penetrance based on age
    if (variant.isAgeSpecific) {
        for (const auto &ageEntry : variant.agePenetrance) {
            if (age >= ageEntry.first) {
                basePenetrance = ageEntry.second;
            }
        }
    } else {
        // For non-age-specific diseases, use a default value
        basePenetrance = 0.5;
    }

    // Adjust penetrance based on genotype
    if (genotype == "0/1") {
        // Heterozygous (carrier)
        basePenetrance *= 0.5; // Example adjustment, could be variant-specific
    } else if (genotype == "1/1") {
        // Homozygous alternate (two copies of the variant)
        basePenetrance *= 1.5; // Example adjustment, could be variant-specific
    } else if (genotype == "0/0") {
        // Homozygous reference (no variant)
        return 0.0; // No risk from this variant
    }

    // Adjust penetrance based on other factors
    for (const auto &factor : userFactors) {
        if (variant.environmentalFactors.find(factor.first) != variant.environmentalFactors.end()) {
            basePenetrance *= variant.environmentalFactors.at(factor.first);
        }
        if (variant.riskFactors.find(factor.first) != variant.riskFactors.end()) {
            basePenetrance *= variant.riskFactors.at(factor.first);
        }
    }

    return basePenetrance;
}

// Function to calculate polygenic risk score (PRS)
float calculatePRS(const vector<Variant> &variants, const vector<Genotype> &genotypes) {
    float prs = 0.0;

    for (const auto &variant : variants) {
        // Find corresponding genotype
        auto it = find_if(genotypes.begin(), genotypes.end(), 
                          [&variant](const Genotype &g) {
                              return g.chromosome == variant.chromosome && g.position == variant.position;
                          });

        if (it != genotypes.end()) {
            // Calculate risk contribution based on genotype and odds ratio
            if (it->genotype == "0/1") {
                prs += log(variant.oddsRatio); // Heterozygous
            } else if (it->genotype == "1/1") {
                prs += 2 * log(variant.oddsRatio); // Homozygous alternate
            }
        }
    }

    // Convert log-odds to probability
    return 1.0 / (1.0 + exp(-prs));
}

bool isValidDNA(const string &dnaSequence) {
    for (char nucleotide : dnaSequence) {
        if (nucleotide != 'A' && nucleotide != 'G' && nucleotide != 'C' && nucleotide != 'T') {
            return false;  // Invalid character found
        }
    }
    return true;  // All characters are valid
}

// Function to detect variants in the user's genome
void detectVariants(const vector<Variant> &variants, const vector<Genotype> &genotypes, int age,
                     const map<string, string> &userFactors) {
    unordered_map<string, float> diseasePenetrance; // Store diseases with their penetrance values
    unordered_map<string, vector<string>> diseaseVariants; // Store variants contributing to each disease

    // Check each variant against the user's genotypes
    for (const auto &variant : variants) {
        auto it = find_if(genotypes.begin(), genotypes.end(), 
                          [&variant](const Genotype &g) {
                              return g.chromosome == variant.chromosome && g.position == variant.position;
                          });

        if (it != genotypes.end() && (it->genotype == "0/1" || it->genotype == "1/1")) {
            // Check quality of the genotype call
            if (it->quality < 20 || it->readDepth < 10) {
                cout << "\n\033[1;33mWarning: Low quality genotype call for " 
                     << variant.chromosome << ":" << variant.position 
                     << " (Quality: " << it->quality << ", Depth: " << it->readDepth 
                     << "). Results may be unreliable.\033[0m" << endl;
            }
            
            // Calculate penetrance for this variant
            float penetrance = calculatePenetrance(variant, age, it->genotype, userFactors);
            
            // Store penetrance and variant info
            if (diseasePenetrance.find(variant.disease) == diseasePenetrance.end()) {
                diseasePenetrance[variant.disease] = penetrance;
            } else {
                // For multiple variants contributing to the same disease, use a more sophisticated model
                // For example, combine penetrance values based on the disease model (e.g., additive, multiplicative)
                diseasePenetrance[variant.disease] = 1.0 - (1.0 - diseasePenetrance[variant.disease]) * (1.0 - penetrance);
            }
            
            // Store variant info
            diseaseVariants[variant.disease].push_back(variant.chromosome + ":" + to_string(variant.position) +
                                                      " " + variant.ref + ">" + variant.alt + " " + it->genotype);
        }
    }

    // Calculate polygenic risk scores for common diseases
    unordered_map<string, float> diseasePRS;
    vector<string> commonDiseases = {"Type 2 Diabetes", "Coronary Artery Disease", "Hypertension"}; // Example common diseases
    
    for (const auto &disease : commonDiseases) {
        // Filter variants for this disease
        vector<Variant> diseaseVariants;
        for (const auto &variant : variants) {
            if (variant.disease == disease) {
                diseaseVariants.push_back(variant);
            }
        }
        
        // Calculate PRS
        if (!diseaseVariants.empty()) {
            diseasePRS[disease] = calculatePRS(diseaseVariants, genotypes);
        }
    }

    // Write data for plotting
    ofstream plotData("plot_data.txt");
    if (!plotData.is_open()) {
        cerr << "Error: Could not open file to write plot data." << endl;
        return;
    }

    // Write headers
    plotData << "Disease\tPenetrance\tLowerCI\tUpperCI\n";

    for (const auto &entry : diseasePenetrance) {
        // Write disease name and penetrance value to the file, enclosing disease name in quotes
        float penetrance = entry.second;
        float lowerCI = max(0.0f, penetrance - 0.1f); // Example confidence interval
        float upperCI = min(1.0f, penetrance + 0.1f); // Example confidence interval
        
        plotData << "\"" << entry.first << "\"\t" << penetrance * 100 << "\t" << lowerCI * 100 << "\t" << upperCI * 100 << endl;
    }

    // Add PRS scores to the plot
    for (const auto &entry : diseasePRS) {
        float prs = entry.second;
        float lowerCI = max(0.0f, prs - 0.1f); // Example confidence interval
        float upperCI = min(1.0f, prs + 0.1f); // Example confidence interval
        
        plotData << "\"" << entry.first << " (PRS)\"\t" << prs * 100 << "\t" << lowerCI * 100 << "\t" << upperCI * 100 << endl;
    }

    plotData.close();

    // Run Gnuplot and open the image
    system("gnuplot -e \"set terminal png size 2000,800; "
           "set output 'risk_plot.png'; "
           "set title 'Disease Risk Assessment'; "
           "set xlabel 'Disease'; "
           "set ylabel 'Risk (%)'; "
           "set grid; "
           "set style data histogram; "
           "set style histogram errorbars gap 2 lw 1; "
           "set style fill solid border -1; "
           "set xtics rotate by -45; "
           "set rmargin 10; "
           "set lmargin 10; "
           "plot 'plot_data.txt' using 2:3:4:xtic(1) title 'Risk with Confidence Intervals'\"");
    system("start risk_plot.png");

    // Print the results in a clean format
    if (diseasePenetrance.empty() && diseasePRS.empty()) {
        cout << "\n\033[1;33mNo known disease risks were found in the genomic data.\033[0m" << endl;  // Yellow text
    } else {
        cout << "\n\033[1;34mDetected Disease Risks:\033[0m" << endl;  // Blue text
        cout << string(100, '-') << endl;
        cout << left << setw(40) << "Disease" 
             << " | " << left << setw(15) << "Risk (%)" 
             << " | " << left << setw(40) << "Variants" << endl;
        cout << string(100, '-') << endl;

        string highestRiskDisease;
        float highestRisk = 0.0;

        // Print variant-based risks
        for (const auto &entry : diseasePenetrance) {
            cout << left << setw(40) << entry.first << " | ";
            cout << fixed << setprecision(2) << entry.second * 100 << "%" << setw(10) << " | ";
            
            // Print variants contributing to this disease
            for (const auto &variant : diseaseVariants[entry.first]) {
                cout << variant << " ";
            }
            cout << endl;
            
            // Check for the highest risk disease
            if (entry.second > highestRisk) {
                highestRisk = entry.second;
                highestRiskDisease = entry.first;
            }
        }
        
        // Print PRS-based risks
        for (const auto &entry : diseasePRS) {
            cout << left << setw(40) << entry.first + " (PRS)" << " | ";
            cout << fixed << setprecision(2) << entry.second * 100 << "%" << setw(10) << " | ";
            cout << "Polygenic risk score based on multiple variants" << endl;
            
            // Check for the highest risk disease
            if (entry.second > highestRisk) {
                highestRisk = entry.second;
                highestRiskDisease = entry.first + " (PRS)";
            }
        }
        
        cout << "\n\033[1;35mThe disease with the highest risk is:\033[0m " 
             << highestRiskDisease << " with " << fixed << setprecision(2) << highestRisk * 100 << "% risk." << endl;
    }
}

// Main function
int main() {
    // User instructions
    cout << "\033[1;32mWelcome to the Enhanced Genetic Risk Analyzer!\033[0m" << endl; // Green text

    char runAgain;  // Variable to check if the user wants to run the program again

    do {
        cout << "\nThis tool analyzes your genetic data to assess potential health risks.\n" << endl;

        // Get input file type
        cout << "\033[1;36mSelect input file type:\033[0m\n";
        cout << "1. VCF file (Variant Call Format)\n";
        cout << "2. FASTA file\n";
        cout << "3. Manual DNA sequence input\n";
        cout << "Enter your choice (1-3): ";
        
        int fileChoice;
        while (!(cin >> fileChoice) || fileChoice < 1 || fileChoice > 3) {
            cout << "\033[1;31mInvalid input. Please enter a number between 1 and 3:\033[0m ";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        
        vector<Genotype> genotypes;
        string dnaSequence;
        
        if (fileChoice == 1) {
            // VCF file input
            string vcfFile;
            cout << "\033[1;36mEnter the path to the VCF file:\033[0m ";
            cin >> vcfFile;
            
            genotypes = parseVCF(vcfFile);
            if (genotypes.empty()) {
                cout << "\033[1;31mError: No valid genotypes found in the VCF file. Please check the file and try again.\033[0m" << endl;
                continue;
            }
        } else if (fileChoice == 2) {
            // FASTA file input
            string fastaFile;
            cout << "\033[1;36mEnter the path to the FASTA file:\033[0m ";
            cin >> fastaFile;
            
            dnaSequence = parseFASTA(fastaFile);
            if (dnaSequence.empty()) {
                cout << "\033[1;31mError: No valid sequence found in the FASTA file. Please check the file and try again.\033[0m" << endl;
                continue;
            }
            
            // For FASTA input, we need to create genotypes by matching with known variants
            // This is a simplified approach, in reality, alignment and variant calling would be more complex
            vector<Variant> variants = parseJSON("disease_data.json");
            for (const auto &variant : variants) {
                // Check if the variant's reference allele is in the DNA sequence
                size_t pos = dnaSequence.find(variant.ref);
                if (pos != string::npos && pos <= variant.position && pos + variant.ref.length() > variant.position) {
                    // Check if the alternate allele is different from the reference
                    if (variant.alt != dnaSequence.substr(variant.position - pos, variant.alt.length())) {
                        genotypes.push_back({variant.chromosome, variant.position, variant.ref, variant.alt, "0/1", 30.0, 20});
                    } else {
                        genotypes.push_back({variant.chromosome, variant.position, variant.ref, variant.alt, "0/0", 30.0, 20});
                    }
                }
            }
        } else {
            // Manual DNA sequence input
            cout << "\033[1;36mEnter the DNA sequence (only A, G, C, T):\033[0m ";
            cin >> dnaSequence;

            // Validate DNA sequence
            while (!isValidDNA(dnaSequence)) {
                cout << "\033[1;31mInvalid DNA sequence. Please enter only A, G, C, T:\033[0m ";
                cin >> dnaSequence;
            }
            
            // For manual input, we need to create genotypes by matching with known variants
            // This is a simplified approach, in reality, alignment and variant calling would be more complex
            vector<Variant> variants = parseJSON("disease_data.json");
            for (const auto &variant : variants) {
                // Position-based matching instead of simple string search
                if (variant.position < dnaSequence.length()) {
                    string userAllele = dnaSequence.substr(variant.position, variant.ref.length());
                    if (userAllele == variant.ref) {
                        // User has reference allele
                        genotypes.push_back({variant.chromosome, variant.position, variant.ref, variant.alt, "0/0", 30.0, 20});
                    } else if (userAllele == variant.alt) {
                        // User has alternate allele
                        genotypes.push_back({variant.chromosome, variant.position, variant.ref, variant.alt, "1/1", 30.0, 20});
                    } else {
                        // Potentially heterozygous or complex variant
                        genotypes.push_back({variant.chromosome, variant.position, variant.ref, variant.alt, "0/1", 30.0, 20});
                    }
                }
            }
        }

        int age;
        cout << "\033[1;36mEnter your age:\033[0m ";
        while (!(cin >> age) || (age <= 0 || age > 100)) {
            cout << "\033[1;31mInvalid input. Please enter a valid age (1-100):\033[0m ";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        
        // Get additional user factors
        map<string, string> userFactors;
        
        cout << "\033[1;36mDo you have a family history of any diseases? (y/n):\033[0m ";
        char familyHistory;
        cin >> familyHistory;
        if (familyHistory == 'y' || familyHistory == 'Y') {
            cout << "\033[1;36mEnter the disease name:\033[0m ";
            string disease;
            cin.ignore(); // Clear newline from previous input
            getline(cin, disease);
            userFactors["FamilyHistory"] = disease;
        }
        
        cout << "\033[1;36mDo you smoke? (y/n):\033[0m ";
        char smoking;
        cin >> smoking;
        userFactors["Smoking"] = (smoking == 'y' || smoking == 'Y') ? "Yes" : "No";
        
        // Load variants from JSON file
        vector<Variant> variants = parseJSON("disease_data.json");

        // Perform variant detection and risk assessment
        detectVariants(variants, genotypes, age, userFactors);

        cout << "\n\033[1;33mDo you want to run the program again? (y/n): \033[0m";
        cin >> runAgain;
        
        // Input validation for run again
        while (runAgain != 'y' && runAgain != 'n') {
            cout << "\033[1;31mInvalid input. Please enter 'y' or 'n': \033[0m";
            cin >> runAgain;
        }

    } while (runAgain == 'y'); // Repeat if user enters 'y'

    cout << "\033[1;32mThank you for using the Enhanced Genetic Risk Analyzer! Goodbye.\033[0m" << endl; // Goodbye message
    return 0;
}
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <tuple>
#include <string>

class pid_gate {
public:
    // Default constructor
    pid_gate() = default;

    // Constructor with parameters
    pid_gate(int frsz, int frsa, int fragz, int fraga, float frszcentre, float frszsigma, float frsaqcentre, float frsaqsigma, float fragzcentre, float fragzsigma, float fragaqcentre, float fragaqsigma)
        : frsz(frsz), frsa(frsa), fragz(fragz), fraga(fraga), frszcentre(frszcentre), frszsigma(frszsigma), frsaqcentre(frsaqcentre), frsaqsigma(frsaqsigma), fragzcentre(fragzcentre), fragzsigma(fragzsigma), fragaqcentre(fragaqcentre), fragaqsigma(fragaqsigma) {}

    // Custom key type for std::map
    using KeyType = std::tuple<int, int, int, int>;

    // Constructor with KeyType
    pid_gate(const KeyType& key)
        : frsz(std::get<0>(key)), frsa(std::get<1>(key)), fragz(std::get<2>(key)), fraga(std::get<3>(key)), frszcentre(0.0), frszsigma(0.0), frsaqcentre(0.0), frsaqsigma(0.0), fragzcentre(0.0), fragzsigma(0.0), fragaqcentre(0.0), fragaqsigma(0.0) {}

    // Accessors for respective protected variables
    int getFrsz() const { return frsz; }
    int getFrsa() const { return frsa; }
    int getFragz() const { return fragz; }
    int getFraga() const { return fraga; }
    float getFrszCentre() const { return frszcentre; }
    float getFrszSigma() const { return frszsigma; }
    float getFrsaqCentre() const { return frsaqcentre; }
    float getFrsaqSigma() const { return frsaqsigma; }
    float getFragzCentre() const { return fragzcentre; }
    float getFragzSigma() const { return fragzsigma; }
    float getFragaqCentre() const { return fragaqcentre; }
    float getFragaqSigma() const { return fragaqsigma; }

    // Function to generate the key for the std::map
    KeyType getKey() const {
        return std::make_tuple(frsz, frsa, fragz, fraga);
    }

    // Friend functions to read and write CSV files
    friend std::vector<pid_gate> readCSV(const std::string& filename);
    friend void writeCSV(const std::string& filename, const std::vector<pid_gate>& pidData);

private:
    int frsz;
    int frsa;
    int fragz;
    int fraga;
    float frszcentre;
    float frszsigma;
    float frsaqcentre;
    float frsaqsigma;
    float fragzcentre;
    float fragzsigma;
    float fragaqcentre;
    float fragaqsigma;
};
// Function to write a CSV file from a vector of pid objects
//void writeCSV(const std::string& filename, const std::vector<pid_gate>& pidData) {
void writeCSV(const TString& filename, const std::vector<pid_gate>& pidData) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Write the column labels in the first line
    file << "FRSZ,FRSA,FragZ,FragA,FRSZcentre,FRSZsigma,FRSAQcentre,FRSAQsigma,FragZcentre,FragZsigma,FragAQcentre,FragAQsigma\n";

    // Write the data for each pid object
    for (const auto& p : pidData) {
        file << p.getFrsz() << "," << p.getFrsa() << "," << p.getFragz() << "," << p.getFraga()
             << "," << p.getFrszCentre() << "," << p.getFrszSigma() << "," << p.getFrsaqCentre() << "," << p.getFrsaqSigma()
             << "," << p.getFragzCentre() << "," << p.getFragzSigma() << "," << p.getFragaqCentre() << "," << p.getFragaqSigma() << "\n";
    }

    file.close();
}

// Function to read a CSV file and store data in a vector of pid objects
std::vector<pid_gate> readCSV(const TString& filename) {
    std::vector<pid_gate> pidData;
    std::map<pid_gate::KeyType, pid_gate> pidMap;

    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return pidData;
    }

    // Skip the header line
    std::string header;
    std::getline(file, header);

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int frsz, frsa, fragz, fraga;
        float frszcentre, frszsigma, frsaqcentre, frsaqsigma, fragzcentre, fragzsigma, fragaqcentre, fragaqsigma;

        // Parse the data from the line
        char comma; // To read the commas between values
        if (!(iss >> frsz >> comma >> frsa >> comma >> fragz >> comma >> fraga >> comma >> frszcentre >> comma >> frszsigma >> comma >> frsaqcentre >> comma >> frsaqsigma >> comma >> fragzcentre >> comma >> fragzsigma >> comma >> fragaqcentre >> comma >> fragaqsigma)) {
            std::cerr << "Error: Invalid data format in line: " << line << std::endl;
            continue;
        }

        pid_gate p(frsz, frsa, fragz, fraga, frszcentre, frszsigma, frsaqcentre, frsaqsigma, fragzcentre, fragzsigma, fragaqcentre, fragaqsigma);
        pidMap[p.getKey()] = p; // Store the object in the map
    }

    file.close();

    // Convert the std::map to std::vector
    for (const auto& kv : pidMap) {
        pidData.push_back(kv.second);
    }

    return pidData;
}


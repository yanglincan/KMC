#ifndef __KMCDOC_H__
#define __KMCDOC_H__
#include "kmc.h"
#include <fstream>

//#define KMCDOC_API __declspec(dllexport)

using namespace std;

class KmcDoc
{
public:
    KmcDoc(const string& strSimuFile);
    ~KmcDoc();
    void ScanFileAndDoSimu();
    void Exit(const string& strErrorTip);
    vector<Reaction*>& GetReactionSet() {
        return *m_pReactionSet;
    }
    const double& GetVolume() const { 
        return m_Volume; 
    }
    vector<SubstanceItem*>* GetSubtanceItemLibrary() { 
        return m_pSubtanceItemLibrary; 
    }
    const string& GetTerminatrName() const {
        return m_TerminatrName; 
    }
    const double& GetTerminatrCove() const {
        return m_TerminationConv; 
    }
    vector<StepsAdder*>& GetStepsAdder() { 
        return *m_pStepsAddArray; 
    }
    const double& GetIntervalTime() const {
        return m_IntervalTime; 
    }
    const string& GetSimuFileName() const {
        return m_strSimuFile;
    }
    void SetSubtanceRole(const string& name, const SubtanceRole& role);
    void SetSubtanceItemLibrary(const string& name, const int& number, const double& weight);
    void SetReactionOrder(const int& mark, const OrderController& order);
    void SetVolume(const double& volue) {
        m_Volume = volue; 
    }
    void SetIntervalTime(const double& time) { 
        m_IntervalTime = time; 
    }
    void SetTermination(const double& terminationconv) {
        m_TerminationConv = terminationconv; 
    }
    void SetTerminatrName(const string& strName) { 
        m_TerminatrName = strName; 
    }
    void SetMultiFeed(StepsAdder* const pAdder) { 
        m_pStepsAddArray->push_back(pAdder);
    }
private:
    void DeleteReactionSetList();
    void DeleteSubtanceItemLibraryVector();
    void DeleteStepsAddArrayVector();
    void DeleteGillespieandAnalysisTool();
    bool InitSubtanceItemLibrary();
    bool CheckParaValuable();
    bool CheckElementRole();
    bool CheckAcountofReactant();
    bool CheckTermination();
    void PrintAllPara();
    void UpdateReactionSetInfo();
private:
    void ScanFileSection(const string& strBegin, const string& strEnd, void(*Analysis)(const string& strLine, KmcDoc* pCurrent));
    friend void ScanReactionSet(const string& strLine, KmcDoc* pCurrent);
    friend void ScanNumAndMw(const string& strLine, KmcDoc* pCurrent);
    friend void ScanSubtanceRole(const string& strLine, KmcDoc* pCurrent);
    friend void ScanReactionOrder(const string& strLine, KmcDoc* pCurrent);
    friend void ScanKMCPara(const string& strLine, KmcDoc* pCurrent);
    friend void ScanMultiFeed(const string& strLine, KmcDoc* pCurrent);
public:
    void SplitString(const string& s, vector<string>& v, const string& c);
private:
    ifstream m_SimuFile;
    string m_strSimuFile;
    Gillespie* m_pGillespie;
    AnalysisTool* m_pAnalysisTool;
    vector<Reaction*>* m_pReactionSet;     
    vector<SubstanceItem*>* m_pSubtanceItemLibrary;   
    vector<StepsAdder*>* m_pStepsAddArray;	
    double m_Volume;
    double m_TerminationConv;
    string m_TerminatrName;	
    double m_IntervalTime;
};
#endif
#include "kmcdoc.h"
#include "reaction.h"
#include  "gillespie.h"
#include "analysistool.h"
using namespace std;

void ScanReactionSet(const string& strLine, KmcDoc* pCurrent)
{
    Reaction* pReaction = new Reaction(strLine);
    false == pReaction->IsReactionSucceeded() ? delete pReaction : pCurrent->GetReactionSet().push_back(pReaction);
}

void KmcDoc::SplitString(const string& s, vector<string>& v, const string& c)
{
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2 - pos1));
        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

void ScanSubtanceRole(const string& strLine, KmcDoc* pCurrent)
{
    vector<string> vectorStr;
    pCurrent->SplitString(strLine, vectorStr, " ");
    size_t iPos_LEFT_Bracket = vectorStr[0].find('[');
    size_t iPos_RIGHT_Bracket = vectorStr[0].find(']');
    size_t count = iPos_RIGHT_Bracket - iPos_LEFT_Bracket - 1;
    pCurrent->SetSubtanceRole(vectorStr[0].substr(iPos_LEFT_Bracket + 1, count),
        (SubtanceRole)atoi(vectorStr[1].data()));
    vectorStr.clear();
}

void ScanNumAndMw(const string& strLine, KmcDoc* pCurrent)
{
    vector<string> vectorStr;
    pCurrent->SplitString(strLine, vectorStr, " ");
    size_t iPos_LEFT_Bracket = vectorStr[0].find('[');
    size_t iPos_RIGHT_Bracket = vectorStr[0].find(']');
    size_t count = iPos_RIGHT_Bracket - iPos_LEFT_Bracket - 1;
    pCurrent->SetSubtanceItemLibrary(vectorStr[0].substr(iPos_LEFT_Bracket + 1, count),
        atoi(vectorStr[1].data()), atof(vectorStr[2].data()));
    vectorStr.clear();
}

void ScanReactionOrder(const string& strLine, KmcDoc* pCurrent)
{
    vector<string> vectorStr;
    pCurrent->SplitString(strLine, vectorStr, " ");
    size_t iPos_MARK_Bracket = vectorStr[0].find('#');
    if (string::npos == iPos_MARK_Bracket)
        return;
    pCurrent->SetReactionOrder(atoi(vectorStr[0].substr(iPos_MARK_Bracket + 1, 2).data()),
        OrderController(atof(vectorStr[1].data()), atof(vectorStr[2].data())));
    vectorStr.clear();
}

void ScanKMCPara(const string& strLine, KmcDoc* pCurrent)
{
    vector<string> vectorStr;
    pCurrent->SplitString(strLine, vectorStr, " ");
    if ("Volume" == vectorStr[0])
        pCurrent->SetVolume(atof(vectorStr[1].data()));
    if ("Sampling_Time" == vectorStr[0])
        pCurrent->SetIntervalTime(atof(vectorStr[1].data()));
    if ("Termination_Condition" == vectorStr[0])
    {
        vector<string> vectorStr2;
        pCurrent->SplitString(vectorStr[1], vectorStr2, ">");
        size_t iPos_LEFT_Bracket = vectorStr2[0].find('[');
        size_t iPos_RIGHT_Bracket = vectorStr2[0].find(']');
        size_t count = iPos_RIGHT_Bracket - iPos_LEFT_Bracket - 1;
        pCurrent->SetTerminatrName(vectorStr2[0].substr(iPos_LEFT_Bracket + 1, count));
        pCurrent->SetTermination(atof(vectorStr2[1].data()));
        vectorStr2.clear();
    }
    vectorStr.clear();
}

void ScanMultiFeed(const string& strLine, KmcDoc* pCurrent)
{
    vector<string> vectorStr;
    pCurrent->SplitString(strLine, vectorStr, " ");
    size_t iPos_MARK_Bracket = vectorStr[0].find('#');
    if (string::npos == iPos_MARK_Bracket)
        return;
    size_t iPos_LEFT_Bracket = vectorStr[2].find('[');
    size_t iPos_RIGHT_Bracket = vectorStr[2].find(']');
    size_t count = iPos_RIGHT_Bracket - iPos_LEFT_Bracket - 1;
    pCurrent->SetMultiFeed(new  StepsAdder(
        vectorStr[2].substr(iPos_LEFT_Bracket + 1, count),
        atoi(vectorStr[0].substr(iPos_MARK_Bracket + 1, 2).data()),
        atof(vectorStr[1].data()),
        atoi(vectorStr[3].data())));
    vectorStr.clear();
}

KmcDoc::KmcDoc(const string& strSimuFile):
    m_strSimuFile(strSimuFile),
    m_pGillespie(NULL),
    m_pAnalysisTool(NULL),
    m_pReactionSet(NULL),
    m_pSubtanceItemLibrary(NULL),
    m_pStepsAddArray(NULL),
    m_Volume(0), 
    m_TerminationConv(0),
    m_TerminatrName(""),
    m_IntervalTime(0)
{
    m_pReactionSet = new vector<Reaction*>;
    m_pSubtanceItemLibrary = new  vector<SubstanceItem*>;
    m_pStepsAddArray = new vector<StepsAdder*>;
    m_pGillespie = new Gillespie(this);
    m_pAnalysisTool = new AnalysisTool;
    m_SimuFile.open(m_strSimuFile, ios::in);
    if (!m_SimuFile.is_open())
        Exit("Fail to open file!");
    ScanFileAndDoSimu();
}

void KmcDoc::DeleteReactionSetList()
{
    if (NULL == m_pReactionSet)
        return;
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
        delete *itor;
    delete m_pReactionSet;
    m_pReactionSet = NULL;
}
void KmcDoc::DeleteSubtanceItemLibraryVector()
{
    if (NULL == m_pSubtanceItemLibrary)
        return;
    for (auto itor = m_pSubtanceItemLibrary->begin(); itor != m_pSubtanceItemLibrary->end(); itor++)
        delete *itor;
    delete m_pSubtanceItemLibrary;
    m_pSubtanceItemLibrary = NULL;
}
void KmcDoc::DeleteStepsAddArrayVector()
{
    if (NULL == m_pStepsAddArray)
        return;
    for (auto itor = m_pStepsAddArray->begin(); itor != m_pStepsAddArray->end(); itor++)
        delete *itor;
    delete m_pStepsAddArray;
    m_pStepsAddArray = NULL;
}
void KmcDoc::DeleteGillespieandAnalysisTool()
{
    if (NULL != m_pGillespie)
    {
        delete m_pGillespie;
        m_pGillespie = NULL;
    }
    if (NULL != m_pAnalysisTool)
    {
        delete m_pAnalysisTool;
        m_pAnalysisTool = NULL;
    }
}
KmcDoc::~KmcDoc()
{
    m_SimuFile.close();
    DeleteReactionSetList();
    DeleteSubtanceItemLibraryVector();
    DeleteStepsAddArrayVector();
    DeleteGillespieandAnalysisTool();
}

void KmcDoc::SetSubtanceRole(const string& name, const SubtanceRole& role)
{
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
    {
        for (auto itor2 = (*itor)->GetElementVector().begin(); itor2 != (*itor)->GetElementVector().end(); itor2++)
        {
            if (name == (*itor2)->elementName)
                (*itor2)->elementRole = role;
        }
    }
}

bool KmcDoc::InitSubtanceItemLibrary()
{
    if (0 == m_pReactionSet->size())
        return false;
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
    {
        for (auto itor2 = (*itor)->GetElementVector().begin(); itor2 != (*itor)->GetElementVector().end(); itor2++)
        {
            bool nFindFlag = false;
            for (auto itor3 = m_pSubtanceItemLibrary->begin(); itor3 != m_pSubtanceItemLibrary->end(); itor3++)
            {
                if ((*itor2)->elementName == (*itor3)->substanceName)
                {
                    nFindFlag = true;
                    break;
                }
                else
                    nFindFlag = false;
            }
            if (!nFindFlag)
                m_pSubtanceItemLibrary->push_back(new SubstanceItem((*itor2)->elementName));
        }
    }
    return true;
}

void KmcDoc::SetSubtanceItemLibrary(const string& name, const int& number, const double& weight)
{
    for (auto itor = m_pSubtanceItemLibrary->begin(); itor != m_pSubtanceItemLibrary->end(); itor++)
    {
        if (name == (*itor)->substanceName)
        {
            (*itor)->substanceAmount = number; 
            (*itor)->substanceMolecularWeight= weight;
        }
    }
}

void KmcDoc::SetReactionOrder(const int& mark, const OrderController& order)
{
    //反应级数赋值
    int iMark = 0;
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
    {
        if (iMark == (mark - 1))
            (*itor)->SetRelationBetween_kMC_k(order);
        ++iMark;
    }
}

bool KmcDoc::CheckElementRole()
{
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
    {
        for (auto itor2 = (*itor)->GetElementVector().begin(); itor2 != (*itor)->GetElementVector().end(); itor2++)
        {
            if (UNSPECIFIED == (*itor2)->elementRole)
                return false;
        }
    }
    return true;
}

bool KmcDoc::CheckAcountofReactant()
{
    for (auto itor = m_pSubtanceItemLibrary->begin(); itor != m_pSubtanceItemLibrary->end(); itor++)
    {
        for (auto itor2 = m_pReactionSet->begin(); itor2 != m_pReactionSet->end(); itor2++)
        {
            for (auto itor3 = (*itor2)->GetElementVector().begin(); itor3 != (*itor2)->GetElementVector().end(); itor3++)
            {
                if ((*itor)->substanceName == (*itor3)->elementName)
                {
                    if (REACTANT == (*itor3)->elementType)
                        (*itor)->substanceisReactant = true;
                    if (PRODUCT == (*itor3)->elementType)
                        (*itor)->substanceisProduct = true;
                }
            }
        }
        if ((*itor)->substanceisReactant && !(*itor)->substanceisProduct && !(*itor)->substanceAmount)
            return false;
    }
    return true;
}

bool KmcDoc::CheckTermination()
{
    if (m_TerminationConv > 1 || m_TerminationConv < 0)
        return false;
    int nFlag = 0;
    for (auto itor = m_pSubtanceItemLibrary->begin(); itor != m_pSubtanceItemLibrary->end(); itor++)
    {
        if ((*itor)->substanceName == m_TerminatrName)
        {
            nFlag = 1;
            break;
        }
    }
    if (0 == nFlag)
        return false;
    return true;
}

bool KmcDoc::CheckParaValuable()
{
    if (false == CheckElementRole())
        return false;
    if (false == CheckAcountofReactant())
        return false;
    if (false == CheckTermination())
        return false;
    if (0 == m_Volume || 0 > m_Volume)
        return false;
    if (m_IntervalTime < 0)
        return false;
    return true;
}

void KmcDoc::PrintAllPara()
{
    cout << endl;
    cout << "{Elementary reaction}:" << endl;
    for (auto itor = m_pReactionSet->begin(); itor != m_pReactionSet->end(); itor++)
         cout << "\t" <<**itor;
    cout << endl << "{Molecular Property}:" << endl;
    for (auto itor = m_pSubtanceItemLibrary->begin(); itor != m_pSubtanceItemLibrary->end(); itor++)
        cout << "\t" << **itor;
    cout << endl << "{Multistep Feeding}:" << endl;
    for (auto itor = m_pStepsAddArray->begin(); itor != m_pStepsAddArray->end(); itor++)
        cout << "\t" << **itor;
    cout << endl << "{Other parameters}:" << endl;
    cout << "\t" << "Simulation Volume: " << m_Volume << endl;
    cout << "\t" << "Sampling time: " << m_IntervalTime <<" s"<< endl;
    cout << "\t" << "Termination condition: " << m_TerminatrName << " > " << m_TerminationConv <<endl;
    cout << endl;
}

void KmcDoc::Exit(const string& strErrorTip)
{
    cout << endl << strErrorTip << endl << endl << "Input \"exit\" to exit!" << endl;
    string str("");
    while ("exit" != str)
        cin >> str;
    exit(1);
}

void KmcDoc::ScanFileSection(const string& strBegin, const string& strEnd, void(*Analysis)(const string& strLine, KmcDoc* pCurrent))
{
    char buffer[256] = { '\0' };
    m_SimuFile.seekg(0, std::ios::beg);
    while (!m_SimuFile.eof())
    {
        m_SimuFile.getline(buffer, 100);
        string strLine(buffer);
        if (string::npos != strLine.find(strBegin))
        {
            while (!m_SimuFile.eof())
            {
                m_SimuFile.getline(buffer, 100);
                strLine = buffer;
                if (string::npos != strLine.find(strEnd))
                    return;
                Analysis(strLine, this);
            }
        }
    }
}

void KmcDoc::UpdateReactionSetInfo()
{
    for (auto itor_r = m_pReactionSet->begin(); itor_r != m_pReactionSet->end(); itor_r++)
    {
        string strAddUnit = "", strEndUnit = "", strInitiator = "", strInitiatSpecie = "";
        for (auto itor = (*itor_r)->GetElementVector().begin(); itor != (*itor_r)->GetElementVector().end(); itor++)
        {
            if (REACTANT == (*itor)->elementType)
            {
                if (INITIATOR == (*itor)->elementRole)
                    strInitiator = (*itor)->elementName;
                else if (MONOMER == (*itor)->elementRole)
                    strAddUnit = (*itor)->elementName;
                else if (SPECIE == (*itor)->elementRole)
                    strInitiatSpecie = (*itor)->elementName;
            }
            else if(PRODUCT == (*itor)->elementType)
            {
                if (SPECIE == (*itor)->elementRole)
                    strEndUnit = (*itor)->elementName;
            }
        }
        (*itor_r)->SetReactionInfo(ReactionInfo(strAddUnit, strEndUnit, strInitiator, strInitiatSpecie));
    }
}

void KmcDoc::ScanFileAndDoSimu()
{
    ScanFileSection("{Elementary Reactions}", "{End Reac}", ScanReactionSet);
    ScanFileSection("{Molecular Property}", "{End Prop}", ScanSubtanceRole);
    if (false == InitSubtanceItemLibrary())
        return;
    UpdateReactionSetInfo();
    ScanFileSection("{Molecular Number and Weight}", "{End Num}", ScanNumAndMw);
    ScanFileSection("{Modify Reaction Order}", "{End Order}", ScanReactionOrder);
    ScanFileSection("{KMC Parameters}", "{End Para}", ScanKMCPara);
    ScanFileSection("{Multistep Feeding}", "{End Multi}", ScanMultiFeed);
    if (false == CheckParaValuable())
        Exit("Error parameter!");
    PrintAllPara();
    cout << "### Confirm the parameters are correct ###" << endl << endl;
    system("pause");
    
    time_t tBegin = time(0);            //计时开始
    m_pGillespie->StartSimulation();    //开始模拟
    m_pAnalysisTool->SetAllChainandKmcDoc(m_pGillespie->GetAllChain(), this);   //将 模拟结果 传入 m_pAnalysisTool
    m_pAnalysisTool->StartAnalysis();   //开始分析
    time_t tEnd = time(0);
    cout << endl << "-> Simulation and analysis totally take " << tEnd - tBegin << " s. <-" << endl << endl;
}
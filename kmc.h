#ifndef __KMC_H__
#define __KMC_H__
#include <string>
#include <vector>
#include <deque>
//#include <list>
#include <stack>
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <algorithm>
//#include <functional>

using namespace std;
enum SubtanceType { REACTANT, PRODUCT };
enum SubtanceRole { UNSPECIFIED, MONOMER, SPECIE, INITIATOR };
enum SyntaxError { SYNTAX_PASS, NON_SMALL_BRACKET , NON_BRACKET, NON_BIG_BRACKET, NON_ARROW, NON_RATE, BRACKET_ERROR,
	NON_NAME, NON_STOICHIOMETRIC_NUMBER, NO_REACTION, MORE_THAN_THREE, ONLY_REACTANT_OR_PRODUCT};
enum ReactionOrder { _UNSPECIFIED_TYPE, _ZERO_TYPE, _A_TYPE, _AB_TYPE, _2A_TYPE, _ABC_TYPE, _A2B_YTPE, _3A_TYPE };

class Reaction;
class KmcDoc;
class Gillespie;
class AnalysisTool;
class KmcException;

struct ReactionInfo {
	string _strAddUnit; 
	string _strEndUnit;
	string _strInitiator;
	string _strInitiatSpecie;
	ReactionInfo(const string add, const  string end, const string initiator, const string specie):
		_strAddUnit(add),
		_strEndUnit(end),
		_strInitiator(initiator),
		_strInitiatSpecie(specie){}
	ReactionInfo() {}
};

struct SpecieVector
{
	string _SpecieName;
	vector<int>* _pVector;
	SpecieVector(string name = "", vector<int>* p = NULL) :_SpecieName(name), _pVector(p) {}
	~SpecieVector() {
		if (NULL != _pVector)
			delete _pVector;
	}
};

class KmcSyntaxException :public exception
{
public:
	char const* what() const
	{
		switch (m_ErrorCode)
		{
		case NO_REACTION: return "\n【反应方程式】无 反应方程式  ！\n\n 示例：\n\n(1)[A] + (1)[B] > (2)[C] {0.128}\n";
		case NON_SMALL_BRACKET: return "\n【反应方程式】无 化学计量系数标识 “( )” ！\n";
		case NON_BRACKET: return "\n【反应方程式】无 反应物or产物名称标识 “[ ]” ！\n";
		case NON_BIG_BRACKET: return "\n【反应方程式】无 反应速率常数标识符 “{ }” ！\n";
		case NON_ARROW:  return "\n【反应方程式】无 反应指向标识 “>” ！\n";
		case NON_RATE: return "\n【反应方程式】反应速率常数 数值 错误！\n";
		case BRACKET_ERROR: return "\n【反应方程式】“( ) [ ] { }”配对顺序 错误！\n";
		case NON_NAME: return "\n【反应方程式】反应物or产物 名称 为空！\n";
		case NON_STOICHIOMETRIC_NUMBER: return "\n【反应方程式】化学计量系数 数值 为空！\n";
		case MORE_THAN_THREE: return "\n【反应方程式】不支持 反应物化学计量系数之和 > 3 的反应！\n";
		case ONLY_REACTANT_OR_PRODUCT: return "\n【反应方程式】不合理的 反应物 or 产物 分配！\n";
		default:
			break;
		}
		return NULL;
	}
	KmcSyntaxException(SyntaxError error = SYNTAX_PASS) :m_ErrorCode(error) {}
private:
	SyntaxError m_ErrorCode;
};
typedef struct ReactionElement
{
	string elementName;		//物质名称
	SubtanceType elementType;	//REACTANT or PRODUCT
	SubtanceRole elementRole;	//UNSPECIFIED, MONOMER, SPECIE or INITIATOR
	int stoichiometricNumber;	//化学计量数
	ReactionElement(string name = "", SubtanceType type = REACTANT, SubtanceRole role = UNSPECIFIED, int number = 0):
		elementName(name),elementType(type),elementRole(role),stoichiometricNumber(number){}
} ReactionElement;

typedef struct OrderController
{
	//(x,y): k_mc = x*k/[(NaV)^y]，修正反应级数
	double x;
	double y;
	OrderController(double xx = 0, double yy = 0) :x(xx), y(yy) {}
} OrderController;

typedef struct StepsAdder
{
	int stepsNumber;	//分步加料的序号，从1开始连续增加
	double stepsConv;	//在当前转化率时，实施分步加料
	string stepsName;	//分布加料的物质名称
	int stepsAmount;	//分步加料的数量
	bool stepsDone;		//分步加料是否完成的标记
	StepsAdder(string name, int number = 0, double conv = 0, int acount = 0, bool done = false) :
		stepsName(name), 
		stepsNumber(number),
		stepsConv(conv),
		stepsAmount(acount),
		stepsDone(done){}
	friend ostream& operator << (ostream& os, const StepsAdder& adder)
	{
		return os << "#" << adder.stepsNumber << ": Conv(" << adder.stepsConv << ")  "
			<<adder.stepsName << "  " << adder.stepsAmount << endl;
	}
} StepsAdder;

typedef struct SubstanceItem
{
	string substanceName;
	int substanceAmount;
	bool substanceisReactant;
	bool substanceisProduct;
	double substanceMolecularWeight;
	SubstanceItem(string name, int amount = 0, bool isReactant = false, bool isProduct = false, double weight = 0) :
		substanceName(name),
		substanceAmount(amount),
		substanceisReactant(isReactant),
		substanceisProduct(isProduct),
		substanceMolecularWeight(weight){}
	friend ostream& operator << (ostream& os, const SubstanceItem& item)
	{
		os << item.substanceName<<"  ";
		if (item.substanceisReactant) os << "  Reaction";
		if (item.substanceisProduct) os << "  Product";
		os << "  Weight:" << item.substanceMolecularWeight << "  Amount:" << item.substanceAmount << endl;
		return os;
	}
} SubstanceItem;

typedef struct BlockItem
{
	string blockName;
	int blockLength;
	int blockPosInChain;
	size_t blockNumberWithSameName;
	BlockItem(string name, int length, int posinchain, size_t numberwithsamename) :
		blockName(name), blockLength(length),blockPosInChain(posinchain),blockNumberWithSameName(numberwithsamename){}
} BlockItem;

typedef struct Distribution
{
	double fValue;
	int iNumbers;
	Distribution(double value, int num) : fValue(value), iNumbers(num) {};
} Distribution;

typedef struct BlockFloat
{
	string blockName;
	double blockLength;
	double blockPosInChain;
	size_t blockNumberWithSameName;
	int valueableChainNumber;
	double lengthStandardDeviation;
	double positionStandardDeviation;
	vector<double>* pLengthArray;
	vector<double>* pPositionArray;
	vector<Distribution>* pLengthDistributionArray;
	vector<Distribution>* pPositionDistributionArray;
	BlockFloat() :
		blockName(""),
		blockLength(0),
		blockPosInChain(0),
		blockNumberWithSameName(0),
		valueableChainNumber(0),
		lengthStandardDeviation(0),
		positionStandardDeviation(0),
		pLengthArray(NULL),
		pPositionArray(NULL),
		pLengthDistributionArray(NULL),
		pPositionDistributionArray(NULL)
	{
	}
} BlockFloat;

typedef vector<BlockItem*>* pBlockArray;
typedef vector<string>* pPerChain;
typedef vector<BlockFloat*>* pBlockFloatArray;
typedef vector<Distribution*>* pDistributionArray;
typedef int number_with_sameName;
typedef size_t max_count;

struct fBlockInfo
{
	string _name;
	max_count _maxcount;
	vector<vector<BlockItem*>*>* _pfBlocks;
	fBlockInfo(string name, max_count count, vector<vector<BlockItem*>*>* pV) :_name(name),_maxcount(count),_pfBlocks(pV){}
};

#endif

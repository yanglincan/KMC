#ifndef __REACTION_H__
#define __REACTION_H__

#include "kmc.h"
using namespace std;
class Reaction
{
public:
	Reaction(string strReaction = "");
	~Reaction();
	Reaction& operator = (const Reaction& reaction);
	friend ostream& operator << (ostream& os, const Reaction& reaction);
	void SetReactionInfo(const ReactionInfo info) {
		m_ReactionInfo = info;
	}
	Reaction(const Reaction& reaction) { 
		this->operator=(reaction); 
	}
	bool IsReactionSucceeded() const { 
		return m_pElementVector ? true : false; 
	}
	const double& GetReactionRateConstant() const {
		return m_ReactionRateConstant; 
	}
	const ReactionOrder& GetReactionOrder() const { 
		return m_ReactionOrder; 
	}
	const OrderController& GetRelationBetween_kMC_k() const { 
		return m_RelationBetween_kMC_k; 
	}
	const vector<ReactionElement*>& GetElementVector() const {
		return *m_pElementVector; 
	}
	void SetRelationBetween_kMC_k(const OrderController order) { 
		m_RelationBetween_kMC_k = order; 
	}
	const ReactionInfo& GetReactionInfo() const {
		return m_ReactionInfo; 
	}
private:
	void DeleteElementVector();
	void ZeroReactant();
	void OneReactant();
	void TwoReactant();
	void ThreeReactant();
	void MoreThanThreeReactant();
	void CheckSyntaxAndInitReac(string strReaction = "");
	void CheckReactionOrder();
	void ThrowExceptionAndDestruction(KmcSyntaxException error) { DeleteElementVector(); throw error;}
private:
	vector<ReactionElement*>* m_pElementVector;
	ReactionOrder m_ReactionOrder;
	OrderController m_RelationBetween_kMC_k; 
	double m_ReactionRateConstant;
private:
	ReactionInfo m_ReactionInfo;
};
#endif

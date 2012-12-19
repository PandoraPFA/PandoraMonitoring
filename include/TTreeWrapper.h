/**
 *  @file   include/TTreeWrapper.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef TTREE_WRAPPER_H
#define TTREE_WRAPPER_H 1

#include <iostream>
#include <map>
#include <vector>

class TTree;
class TBranch;

//------------------------------------------------------------------------------------------------------------------------------------------

namespace pandora_monitoring
{

typedef std::vector<float>  VectorFloat;
typedef std::vector<double> VectorDouble;
typedef std::vector<int>    VectorInt;

class TTreeWrapper 
{
public:
    // classes to throw in case of error
    class TreeInsertError {};
    class BranchInsertError {};
    class TreeNotFoundError {};
    class BranchNotFoundError {};

    TTreeWrapper();

    ~TTreeWrapper();

    template< typename VarType >
    bool Set(const std::string &treeName, const std::string &branchName, VarType value);

    void Fill(const std::string &treeName);
    void Print(const std::string &treeName) const ;
    void Scan(const std::string &treeName) const ;

    TTree*& GetTree(const std::string& treeName);

    void Clear();  ///< clear tree with name 'treeName' in the TTreeWrapper

    template <typename T>
    class Branch
    {
    public:
        class BadType {};

        Branch(TTree *tree, const std::string &branchName);

        ~Branch();

        void Set(T variable);

    private:
        std::string             m_name;                     ///< 
        TTree                  *m_pTree;                    ///< 
        TBranch                *m_pBranch;                  ///< 
        T                       m_variable;                 ///< 
        bool                    m_isVector;                 ///<
    };

    class BranchHandler
    {
    public:
        BranchHandler(TTree *pTree, const std::string &branchName);
        ~BranchHandler();

        bool Set(float value);
        bool Set(double value);
        bool Set(int value);
        bool Set(VectorFloat *ptr);
        bool Set(VectorDouble *ptr);
        bool Set(VectorInt *ptr);

    private:
        // possible types
        enum BranchType
        {
            BRANCH_FLOAT,
            BRANCH_DOUBLE,
            BRANCH_INT,
            BRANCH_VECTOR_FLOAT,
            BRANCH_VECTOR_DOUBLE,
            BRANCH_VECTOR_INT,
            BRANCH_NO_TYPE_DEFINED
        };

        BranchType              m_branchType;               ///<

        Branch<float>          *m_branchFloat;              ///<
        Branch<double>         *m_branchDouble;             ///<
        Branch<int>            *m_branchInt;                ///<

        Branch<VectorFloat*>   *m_branchVectorFloat;        ///<
        Branch<VectorDouble*>  *m_branchVectorDouble;       ///<
        Branch<VectorInt*>     *m_branchVectorInt;          ///<

        TTree*                  m_tree;                     ///<
        std::string             m_branchName;               ///<
    };

private:
    typedef std::map<const std::string, BranchHandler*> BranchMap;
    typedef std::pair<TTree*, BranchMap*> TreeInfo;
    typedef std::map<const std::string, TreeInfo> TreeMap;

    TreeMap::iterator AddTree(const std::string &treeName);

    template< typename VarType >
    BranchMap::iterator AddBranch(const std::string &treeName, const std::string &branchName);

    TreeMap                     m_treeMap;                  ///< holds treenames and ttrees and Branches 
};

} // namespace pandora_monitoring

#endif // #ifndef TTREE_WRAPPER_H

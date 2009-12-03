/**
 *  @file   TTreeWrapper/include/TTreeWrapper.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef TTREE_WRAPPER_H
#define TTREE_WRAPPER_H 1

#include <iostream>
#include <map>

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

    template <typename T> class Branch {
    public:
	class BadType {};

	Branch(TTree* tree, std::string branchName);
	~Branch();

	void Set( T variable );

    private:
	std::string m_name;
	TTree*   m_pTree;
	TBranch* m_pBranch;
		    
	T m_variable;
	bool isVector;
    };


    class BranchHandler {
	// possible types
	enum BranchType {
	    BRANCH_FLOAT,
	    BRANCH_DOUBLE,
	    BRANCH_INT,
	    BRANCH_VECTOR_FLOAT,
	    BRANCH_VECTOR_DOUBLE,
	    BRANCH_VECTOR_INT,
	    BRANCH_NO_TYPE_DEFINED
	};

	BranchType  m_branchType;

	Branch<float>  *m_branchFloat;
	Branch<double> *m_branchDouble;
	Branch<int>    *m_branchInt;

	Branch<VectorFloat*>    *m_branchVectorFloat;
	Branch<VectorDouble*>   *m_branchVectorDouble;
	Branch<VectorInt*>      *m_branchVectorInt;
	    

	TTree*      m_tree;
	std::string m_branchName;

    public:
	BranchHandler( TTree* tree, std::string branchName );
	~BranchHandler();

	bool Set( float  value );
	bool Set( double value );
	bool Set( int    value );
	bool Set( VectorFloat*  ptr );
	bool Set( VectorDouble* ptr );
	bool Set( VectorInt*    ptr );
    };

 public:

    // classes to throw in case of error
    class TreeInsertError   {};
    class BranchInsertError {};
    class TreeNotFoundError {};


    typedef std::map<std::string, BranchHandler*> BranchMap;
    typedef std::pair<TTree*, BranchMap*> TreeInfo;
    typedef std::map<std::string,TreeInfo> TreeMap;

    TTreeWrapper();
    ~TTreeWrapper();

    template< typename VarType >
	bool Set( std::string treeName, std::string branchName, VarType value );

    void Fill ( std::string treeName );
    void Print( std::string treeName );

    TTree* GetTree( std::string treeName ) const;

 private:

    TreeMap::iterator AddTree( std::string treeName );

    template< typename VarType >
	BranchMap::iterator AddBranch( std::string treeName, std::string branchName );


	    
    TreeMap m_treeMap; // holds treenames and ttrees and Branches 
	    
};



/* template< typename VarType > */
/* bool TTreeWrapper::Set( std::string treeName, std::string branchName, VarType value ) */
/* { */
/*     BranchMap::iterator branchIt; */
/*     try */
/*     { */
/*         branchIt = AddBranch<VarType>( treeName, branchName ); */
/*     } */
/*     catch( BranchInsertError& excpt ) */
/*     { */
/*         throw; */
/*     } */
    
/*     return branchIt->second.Set( value ); */
/* } */


/* //------------------------------------------------------------------------------------------------------------------------------------------ */

/* template< typename VarType > */
/* TTreeWrapper::BranchMap::iterator TTreeWrapper::AddBranch( std::string treeName, std::string branchName ) */
/* { */
/*     TreeMap::iterator treeIt = m_treeMap.end(); */
/*     try */
/*     { */
/*         TreeMap::iterator treeIt = AddTree( treeName ); */
/*     } */
/*     catch( TreeInsertError& excpt ) */
/*     { */
/*         throw; */
/*     } */

/*     BranchMap::iterator branchIt = treeIt->second.second.find( branchName ); */
/*     if( branchIt == treeIt->second.second.end() ) */
/*     { */
/*         TTree* tree = treeIt->second.first; */
/*         std::pair<BranchMap::iterator, bool> itInfo = m_treeMap.insert( BranchMap::value_type( branchName, BranchHandler( tree, branchName ) ) ); */
/*         if( !itInfo.second ) */
/*             throw BranchInsertError(); */

/*         return itInfo.first; */
/*     } */
/*     return branchIt; */
/* } */




} // namespace pandora_monitoring

#endif // #ifndef TTREE_WRAPPER_H

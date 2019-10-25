/*******************************************************************************
    PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
    Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
    Email: profasi@thep.lu.se
    Home Page: http://cbbp.thep.lu.se/activities/profasi/
    Version control (git) : https://trac.version.fz-juelich.de/PROFASI

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License
    (see PROFASI/gpl.txt).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
********************************************************************************/

#ifndef CompositePredicate_HH
#define CompositePredicate_HH
#include <string>

namespace prf
{
    //! An extensible composite-capable simple predicate class
    /**
    * A simple predicate is a function object that returns true or false
    * depending on the state of a variable of a certain type. For instance,
    * a predicate can be written, which takes integer variables as argument,
    * and returns true if they are even. The template argument "T" here, is
    * the type of data object the composite predicate is supposed to work on,
    * integers in the above example.
    * \ingroup pdb_handling
    */
    template <class T>

    class SimplePredicate
    {
    public:
        SimplePredicate() : myid("constant") {myconstantstate=false;}

        virtual ~SimplePredicate() {}

        void set_default_state(bool vl) {
            myconstantstate=vl;
            if (vl) myid=std::string("True"); else myid=std::string("False");
        }

        virtual bool operator()(T x) {return myconstantstate;}

        std::string id() const {return myid;}

        inline void set_id(std::string str) {myid=str;}

    private:
        std::string myid;
        bool myconstantstate;
    };

    //! A composite predicate class
    /**
    * \anchor CompositePredicate
    * An arbitrary logical composite of several predicates can be represented
    * as a binary tree, where each node is either a simple predicate, or a
    * predicate whose value is calculated using the results of the evaluation
    * of the daughter nodes and their combination with a logical operator. This
    * class is meant to represent that abstract idea of composition. The exact
    * nature of the predicates does not matter for the tree structure of the
    * composites. Using the helper class SimplePredicate, it is possible to
    * leave the exact way the predicates handle the data entirely out of
    * consideration when compositing. In other words, one can composite two
    * predicates which do entirely unrelated operations on the data.
    * \ingroup pdb_handling
    *
    */
    template <class T>

    class CompositePredicate
    {
        typedef enum {AND,OR,XOR} CompRule;
    private:
        SimplePredicate<T> *unbranched;
        CompositePredicate<T> *branch1,*branch2;
        CompRule myop;
        bool node_negation;
    public:
        CompositePredicate() : unbranched(NULL), branch1(NULL),branch2(NULL),
                myop(AND), node_negation(false) {}

        ~CompositePredicate() {}

        CompositePredicate(SimplePredicate<T> *, bool negated=false);
        CompositePredicate(const CompositePredicate<T> &);
        CompositePredicate(const CompositePredicate<T> &, bool negated);
        CompositePredicate<T> & operator=(const CompositePredicate<T> &);
        CompositePredicate(CompRule somerule,
                           CompositePredicate<T> &pred1,
                           CompositePredicate<T> &pred2);
        inline void inversion(bool invr) {node_negation=invr;}

        std::string id();
        bool operator()(T x);

        CompositePredicate<T> operator!(void) {
            return CompositePredicate<T>(*this,true);
        }

        CompositePredicate<T> operator||(CompositePredicate<T> &pred) {
            return CompositePredicate<T>(OR,*this,pred);
        }

        CompositePredicate<T> operator&&(CompositePredicate<T> &pred) {
            return CompositePredicate<T>(AND,*this,pred);
        }

        CompositePredicate<T> operator^(CompositePredicate<T> &pred) {
            return CompositePredicate<T>(XOR,*this,pred);
        }
    };

    template <class T>
    CompositePredicate<T>::CompositePredicate(SimplePredicate<T> *sp,
            bool negated)
    {
        unbranched=sp;
        branch1=NULL; branch2=NULL;
        myop=OR;
        node_negation=negated;
    }

    template <class T>
    CompositePredicate<T>::CompositePredicate(const CompositePredicate<T> &cp)
    {
        branch1=cp.branch1;
        branch2=cp.branch2;
        unbranched=cp.unbranched;
        myop=cp.myop;
        node_negation=cp.node_negation;
    }

    template <class T>
    CompositePredicate<T> & CompositePredicate<T>::operator=(
        const CompositePredicate<T> &cp)
    {
        if (this!=&cp) {
            branch1=cp.branch1;
            branch2=cp.branch2;
            unbranched=cp.unbranched;
            myop=cp.myop;
            node_negation=cp.node_negation;
        }

        return *this;
    }

    template<class T>
    bool CompositePredicate<T>::operator()(T x)
    {
        //     prf::cout<<"evaluating "<<id()<<"::operator ()\n";
        if (branch1==NULL || branch2==NULL) {
            if (unbranched==NULL) return node_negation;
            else {
                return node_negation xor(*unbranched)(x);
            }
        }

        switch (myop) {
            case XOR: return node_negation xor
                                 ((*branch1)(x) xor(*branch2)(x));
            case OR: return node_negation xor
                                ((*branch1)(x) || (*branch2)(x));
            case AND:
            default: return node_negation xor
                                ((*branch1)(x) && (*branch2)(x));
        };
    }

    template <class T>
    CompositePredicate<T>::CompositePredicate(const CompositePredicate<T> &cp,
            bool negated)
    {
        branch1=cp.branch1;
        branch2=cp.branch2;
        unbranched=cp.unbranched;
        myop=cp.myop;
        node_negation=(negated xor cp.node_negation);
    }

    template <class T>
    CompositePredicate<T>::CompositePredicate(CompRule somerule,
            CompositePredicate &pred1,
            CompositePredicate &pred2)
    {
        unbranched=NULL;
        branch1=&pred1;
        branch2=&pred2;
        myop=somerule;
        node_negation=false;
    }

    template <class T>
    std::string CompositePredicate<T>::id()
    {
        std::string op;

        if (unbranched!=NULL) {
            std::string nds=unbranched->id();
            if (node_negation) nds="not "+nds;
            if (nds=="not False") return "True";
            else if (nds=="not True") return "False";
            else return nds;
        }

        std::string b1=branch1->id();

        std::string b2=branch2->id();

        std::string res="";

        switch (myop) {
        case OR: {
                op=" or ";
                if (b1=="True" or b2=="True") {
                    res="True";
                } else if (b1=="False" and b2=="False") {
                    res="False";
                } else if (b1=="False") {
                    res=b2;
                } else if (b2=="False") {
                    res=b1;
                } else res="("+b1+op+b2+")";
                break;
            }
        case XOR: {
                op=" xor ";
                if ((b1=="True" and b2=="True") or
                    (b1=="False" and b2=="False")) {
                    res="False";
                } else if (b1=="True") res="not "+b2;
                else if (b1=="False") res=b2;
                else if (b2=="True") res="not "+b1;
                else if (b2=="False") res=b1;
                else res="("+b1+op+b2+")";
                break;
            }
            case AND:
        default: {
                op=" and ";
                if (b1=="False" or b2=="False") {
                    res="False";
                } else if (b1=="True" and b2=="True") {
                    res="True";
                } else if (b1=="True") {
                    res=b2;
                } else if (b2=="True") {
                    res=b1;
                } else res="("+b1+op+b2+")";
                break;
            }
        };

        if (node_negation) {
            if (res.substr(0,3)=="not") res=res.substr(3);
            else res="not "+res;
        }

        return res;
    }

}

#endif

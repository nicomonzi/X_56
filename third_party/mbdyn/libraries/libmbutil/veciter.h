/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Iteratore per vettori */


#ifndef VECITER_H
#define VECITER_H

#ifdef USE_MULTITHREAD
#include <signal.h>
#include <atomic>
#endif /* USE_MULTITHREAD */

#include "myassert.h"

/* GetFirst ritorna true ed assegna a TReturn il primo item del vettore.
 * Si assume che il vettore contenga almeno un termine;
 * GetNext ritorna true ed assegna a TReturn il termine successivo se esiste,
 * altrimenti ritorna false e l'assegnamento a TReturn e' unpredictable. */


template <class T>
class Iter {
public:
	virtual ~Iter(void) { NO_OP; };
	virtual bool bGetFirst(T& TReturn) const = 0;
	virtual bool bGetNext(T& TReturn) const = 0;
};

template<class T>
class VecIter : public Iter<T> {
protected:
	mutable T* pStart;
	mutable T* pCount;
	unsigned iSize;

public:
	VecIter(void) : pStart(0), pCount(0), iSize(0) { NO_OP; };
	VecIter(const T* p, unsigned i) : pStart(p), pCount(p), iSize(i)
	{
		ASSERT(pStart != 0);
		ASSERT(iSize > 0);
	};

	virtual ~VecIter(void)
	{
		NO_OP;
	};

	void Init(const T* p, unsigned i)
	{
		ASSERT(p != NULL);
		ASSERT(i > 0);

		pStart = pCount = const_cast<T *>(p);
		iSize = i;
	};

	inline bool bGetFirst(T& TReturn) const
	{
		ASSERT(pStart != NULL);
		ASSERT(iSize > 0);

		if (pStart == NULL) {
			return false;
		}

		pCount = pStart;
		TReturn = *pStart;

		return true;
	};

	inline bool bGetCurr(T& TReturn) const
	{
		ASSERT(pStart != NULL);
		ASSERT(iSize > 0);
		ASSERT(pCount >= pStart);

		if (pCount == pStart + iSize) {
			return false;
		}

		TReturn = *pCount;

		return true;
	};

	inline bool bGetNext(T& TReturn) const
	{
		ASSERT(pStart != NULL);
		ASSERT(iSize > 0);
		ASSERT(pCount >= pStart);

		++pCount;
		if (pCount == pStart + iSize) {
			return false;
		}

		TReturn = *pCount;

		return true;
	};
};

#ifdef USE_MULTITHREAD
/*
 * The user's class must inherit from InUse to be used by the MT_VecIter
 * the user must reset the inuse flag by using SetInUse() before 
 * concurrently iterating over the array; the iterator provides a 
 * helper routine for this; provide it is called only once and not
 * concurrently.
 */
class InUse {
private:
     mutable volatile std::atomic<bool> inuse;

public:
	InUse(void) : inuse(false) { NO_OP; };
	virtual ~InUse(void) { NO_OP; };

	inline bool bIsInUse(void) const
	{
		/*
		 * If inuse is...
		 * 	true:	leave it as is; return false
		 * 	false:	make it true; return true
		 */
                bool locked = false;
                return inuse.compare_exchange_strong(locked, true);
	};
        inline void ReSetInUse() { inuse.exchange(false); };
};

/* #define DEBUG_VECITER */

template<class T>
class MT_VecIter : public VecIter<T> {
protected:
#ifdef DEBUG_VECITER
	mutable unsigned iCount;
#endif /* DEBUG_VECITER */
        unsigned iOffset; // Each thread should process it's own set of elements defined by iOffset and iStep
        unsigned iStep;
        mutable bool bSecondRound; // After a thread has processed it's own set of elements it may start a second round
                                   // in order to process elements from other threads which were not processed so far.     
public:
        MT_VecIter()
             :iOffset(0), iStep(1), bSecondRound(false)
	{
#ifdef DEBUG_VECITER
             iCount = 0;
#endif
	}

	virtual ~MT_VecIter(void)
	{
		NO_OP;
	}

        void Init(const T* pStart_a, unsigned iSize_a, unsigned iOffset_a = 0, unsigned iStep_a = 1) {
                VecIter<T>::Init(pStart_a, iSize_a);
                this->iOffset = iOffset_a;
                this->iStep = iStep_a;
                bSecondRound = false;
        }
	/* NOTE: it must be called only once */
	void ResetAccessData(void)
	{
		ASSERT(this->pStart != nullptr);
		ASSERT(this->iSize > 0);

		for (unsigned i = 0; i < this->iSize; i++) {
			this->pStart[i]->ReSetInUse();
		}
	}

	inline bool bGetFirst(T& TReturn) const
	{
		ASSERT(this->pStart != nullptr);
		ASSERT(this->iSize > 0);

#ifdef DEBUG_VECITER
		iCount = 0;
#endif /* DEBUG_VECITER */

		this->pCount = this->pStart - iStep + iOffset;
                bSecondRound = false;
                
		return bGetNext(TReturn);
	}

	inline bool bGetCurr(T& TReturn) const
	{
		ASSERT(this->pStart != nullptr);
		ASSERT(this->iSize > 0);
		ASSERT(this->pCount >= this->pStart);

		if (this->pCount >= this->pStart + this->iSize) {
			return false;
		}

		TReturn = *this->pCount;
		/* NOTE: of course, by definition it's already in use */

		return true;
	}

	inline bool bGetNext(T& TReturn) const
	{
		ASSERT(this->pStart != nullptr);
		ASSERT(this->iSize > 0);
		ASSERT(this->pCount >= this->pStart - iStep); 
                ASSERT(this->pCount < this->pStart + this->iSize + iStep);
                
                const unsigned iStepCurr = bSecondRound ? 1u : iStep;
                
		for (this->pCount += iStepCurr; 
                     this->pCount < this->pStart + this->iSize; 
                     this->pCount += iStepCurr) {

                        ASSERT(this->pCount >= this->pStart);
                        
			if ((*this->pCount)->bIsInUse()) {
				TReturn = *this->pCount;
#ifdef DEBUG_VECITER
				iCount++;
#endif /* DEBUG_VECITER */
				return true;
			}
		}

                ASSERT(this->pCount >= this->pStart); 
                ASSERT(this->pCount < this->pStart + this->iSize + iStep);

                if (!bSecondRound) {
                     bSecondRound = true;
                     this->pCount = this->pStart - 1;
                     return bGetNext(TReturn);
                }
                
#ifdef DEBUG_VECITER
		silent_cerr("[" << pthread_self() << "]: total=" << iCount
				<< std::endl);
#endif /* DEBUG_VECITER */
                
		return false;
	}
};

#endif /* USE_MULTITHREAD */

#endif /* VECITER_H */

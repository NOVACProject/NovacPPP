/**
* Contains a defintion of a multiplication function object.
*
* @author		\URL[Stefan Kraus]{http://stefan@00kraus.de} @ \URL[IWR, Image Processing Group]{http://klimt.iwr.uni-heidelberg.de}
* @version		1.0 @ 2002/02/01
*/
#if !defined(MULFUNCTION_H_020201)
#define MULFUNCTION_H_020201

#include "ParamFunction.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#pragma warning (push, 3)

namespace MathFit
{
	/**
	* This object represents a multiplication function.
	*
	* @author		\URL[Stefan Kraus]{http://stefan@00kraus.de} @ \URL[IWR, Image Processing Group]{http://klimt.iwr.uni-heidelberg.de}
	* @version		1.0 @ 2002/02/01
	*/
	class CMulFunction : public IParamFunction
	{
	public:
		/**
		* Creates the object and sets the two operand function objects.
		*
		* @param ipfFirst	The first operand function object.
		* @param ipfSec	The second operand function object.
		*/
		CMulFunction(IParamFunction& ipfFirst, IParamFunction& ipfSec) : mFirst(ipfFirst), mSec(ipfSec)
		{
			BuildLinearParameter();
			BuildNonlinearParameter();
		}

		/**
		* Returns the function's value given by the two multiplicative operands.
		*
		* @param fXValue	The X data point.
		*
		* @return	The value of the function at the given data point.
		*/
		virtual TFitData GetValue(TFitData fXValue)
		{
			return mFirst.GetValue(fXValue) * mSec.GetValue(fXValue);
		}

		/**
		* Returns the first derivative of the function at the given data point.
		*
		* @param fXValue	The X data point at which the first derivation is needed.
		*
		* @return	The slope of the function at the given data point.
		*/
		virtual TFitData GetSlope(TFitData fXValue)
		{
			return mFirst.GetSlope(fXValue) * mSec.GetValue(fXValue) + mFirst.GetValue(fXValue) * mSec.GetSlope(fXValue);
		}

		/**
		* Returns the basis function of the specified linear parameter.
		* A basis function is defined as the term by which the linear parameter is multiplied.
		*
		* @param fXValue	The data point at which the basis function should be determined.
		* @param iParamID	The index within the linear parameter vector of the linear parameter.
		* @param bFixedID	If TRUE the given parameter ID is the parameter ID without all fixed parameter.
		*
		* @return	The basis function in regard to the given linear parameter.
		*/
		virtual TFitData GetLinearBasisFunction(TFitData fXValue, int iParamID, bool bFixedID = true)
		{
			if(bFixedID)
			{
				const int iSize = mFirst.GetLinearParameter().GetSize();
				if(iParamID >= iSize)
					return mSec.GetLinearBasisFunction(fXValue, iParamID - iSize, bFixedID);
				else
					return mFirst.GetLinearBasisFunction(fXValue, iParamID, bFixedID);
			}
			else
			{
				const int iSize = mFirst.GetLinearParameterVector().GetAllSize();
				if(iParamID >= iSize)
					return mSec.GetLinearBasisFunction(fXValue, iParamID - iSize, bFixedID);
				else
					return mFirst.GetLinearBasisFunction(fXValue, iParamID, bFixedID);
			}
		}

		/**
		* Sets the new set of linear parameters.
		* The given parameter vector is split up into the pieces needed by the references.
		*
		* @param vParam	The vector containing the new unfixed parameters.
		*
		* @return	TRUE if successful, FALSE otherwise.
		*/
		virtual bool SetNonlinearParameter(CVector& vParam)
		{
			// keep internal copy
			mNonlinearParams.SetParameters(vParam);

			const int iFirstSize = mFirst.GetNonlinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = vParam.SubVector(0, iFirstSize);
				mFirst.SetNonlinearParameter(temp);
			}

			const int iSecSize = mSec.GetNonlinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = vParam.SubVector(iFirstSize, iSecSize);
				mSec.SetNonlinearParameter(temp);
			}

			return true;
		}

		/**
		* Sets the linear parameters of the function object.
		*
		* @param vLinParam	A vector object containing the new unfixed linear parameters.
		*
		* @return TRUE if successful, FALSE if the vector sizes do not match.
		*/
		virtual bool SetLinearParameter(CVector& vLinParam)
		{
			mLinearParams.SetParameters(vLinParam);

			const int iFirstSize = mFirst.GetLinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = vLinParam.SubVector(0, iFirstSize);
				mFirst.SetLinearParameter(temp);
			}

			const int iSecSize = mSec.GetLinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = vLinParam.SubVector(iFirstSize, iSecSize);
				mSec.SetLinearParameter(temp);
			}

			return true;
		}

		/**
		* Sets the covariance matrix of the nonlinear parameters.
		*
		* @param mCovar	The covariance matrix.
		*/
		virtual void SetNonlinearCovarMatrix(CMatrix& mCovar)
		{
			mNonlinearParams.SetCovarMatrix(mCovar);

			// split up vevctor into linear and nonlinear parts.
			const int iFirstSize = mFirst.GetNonlinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = mCovar.SubMatrix(0, 0, iFirstSize, iFirstSize);
				mFirst.SetNonlinearCovarMatrix(temp);
			}

			const int iSecSize = mSec.GetNonlinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = mCovar.SubMatrix(iFirstSize, iFirstSize, iSecSize, iSecSize);
				mSec.SetNonlinearCovarMatrix(temp);
			}
		}

		/**
		* Sets the correlation matrix of the nonlinear parameters.
		*
		* @param mCorrel	The correlation matrix.
		*/
		virtual void SetNonlinearCorrelMatrix(CMatrix& mCorrel)
		{
			mNonlinearParams.SetCorrelMatrix(mCorrel);

			// split up vevctor into linear and nonlinear parts.
			const int iFirstSize = mFirst.GetNonlinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = mCorrel.SubMatrix(0, 0, iFirstSize, iFirstSize);
				mFirst.SetNonlinearCorrelMatrix(temp);
			}

			const int iSecSize = mSec.GetNonlinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = mCorrel.SubMatrix(iFirstSize, iFirstSize, iSecSize, iSecSize);
				mSec.SetNonlinearCorrelMatrix(temp);
			}
		}

		/**
		* Sets the error of the nonlinear parameters.
		*
		* @param vError	The error vector.
		*/
		virtual void SetNonlinearError(CVector& vError)
		{
			mNonlinearParams.SetError(vError);

			// split up vevctor into linear and nonlinear parts.
			const int iFirstSize = mFirst.GetNonlinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = vError.SubVector(0, iFirstSize);
				mFirst.SetNonlinearError(temp);
			}

			const int iSecSize = mSec.GetNonlinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = vError.SubVector(iFirstSize, iSecSize);
				mSec.SetNonlinearError(temp);
			}
		}

		/**
		* Sets the covariance matrix of the linear parameters.
		*
		* @param mCovar	The covariance matrix.
		*/
		virtual void SetLinearCovarMatrix(CMatrix& mCovar)
		{
			mLinearParams.SetCovarMatrix(mCovar);

			// split up vevctor into linear and nonlinear parts.
			const int iFirstSize = mFirst.GetLinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = mCovar.SubMatrix(0, 0, iFirstSize, iFirstSize);
				mFirst.SetLinearCovarMatrix(temp);
			}

			const int iSecSize = mSec.GetLinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = mCovar.SubMatrix(iFirstSize, iFirstSize, iSecSize, iSecSize);
				mSec.SetLinearCovarMatrix(temp);
			}
		}

		/**
		* Sets the correlation matrix of the linear parameters.
		*
		* @param mCorrel	The correlation matrix.
		*/
		virtual void SetLinearCorrelMatrix(CMatrix& mCorrel)
		{
			mLinearParams.SetCorrelMatrix(mCorrel);

			// split up vevctor into linear and Linear parts.
			const int iFirstSize = mFirst.GetLinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = mCorrel.SubMatrix(0, 0, iFirstSize, iFirstSize);
				mFirst.SetLinearCorrelMatrix(temp);
			}

			const int iSecSize = mSec.GetLinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = mCorrel.SubMatrix(iFirstSize, iFirstSize, iSecSize, iSecSize);
				mSec.SetLinearCorrelMatrix(temp);
			}
		}

		/**
		* Sets the error of the linear parameters.
		*
		* @param vError	The error vector.
		*/
		virtual void SetLinearError(CVector& vError)
		{
			mLinearParams.SetError(vError);

			// split up vevctor into linear and nonlinear parts.
			const int iFirstSize = mFirst.GetLinearParameter().GetSize();
			if(iFirstSize > 0)
			{
				auto temp = vError.SubVector(0, iFirstSize);
				mFirst.SetLinearError(temp);
			}

			const int iSecSize = mSec.GetLinearParameter().GetSize();
			if(iSecSize > 0)
			{
				auto temp = vError.SubVector(iFirstSize, iSecSize);
				mSec.SetLinearError(temp);
			}
		}

		/**
		* Resets the linear parameters to default values.
		* The default implementation sets every parameter to zero.
		*/
		virtual void ResetLinearParameter()
		{
			mFirst.ResetLinearParameter();
			mSec.ResetLinearParameter();
			BuildLinearParameter();
		}

		/**
		* Resets the nonlinear parameters to default values.
		* The default implementation sets every parameter to zero.
		*/
		virtual void ResetNonlinearParameter()
		{
			mFirst.ResetNonlinearParameter();
			mSec.ResetNonlinearParameter();
			BuildNonlinearParameter();
		}

		virtual TFitData GetLinearPenalty(TFitData fChiSquare)
		{
			return mFirst.GetLinearPenalty(fChiSquare) + mSec.GetLinearPenalty(fChiSquare);
		}

		virtual TFitData GetNonlinearPenalty(TFitData fChiSquare)
		{
			return mFirst.GetNonlinearPenalty(fChiSquare) + mSec.GetNonlinearPenalty(fChiSquare);
		}

		virtual void SetFitRange(CVector& vFitRange)
		{
			mFitRange.Copy(vFitRange);
			mFirst.SetFitRange(vFitRange);
			mSec.SetFitRange(vFitRange);
		}

	private:
		/**
		* Builds the linear parameter vector from the operands parameters.
		*/
		void BuildLinearParameter()
		{
			// get the original parameters
			CVector vBuf;
			vBuf.Append(mFirst.GetLinearParameter());
			vBuf.Append(mSec.GetLinearParameter());
			mLinearParams.SetSize(vBuf.GetSize());
			mLinearParams.SetParameters(vBuf);
		}

		/**
		* Builds the nonlinear parameter vector from the operands parameters.
		*/
		void BuildNonlinearParameter()
		{
			CVector vBuf;
			vBuf.Append(mFirst.GetNonlinearParameter());
			vBuf.Append(mSec.GetNonlinearParameter());
			mNonlinearParams.SetSize(vBuf.GetSize());
			mNonlinearParams.SetParameters(vBuf);
		}

		/**
		* Hold the first operand.
		*/
		IParamFunction& mFirst;
		/**
		* Hold the second operand.
		*/
		IParamFunction& mSec;
	};
}

#pragma warning (pop)
#endif

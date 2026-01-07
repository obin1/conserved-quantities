!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Matlab Gateway for the Sparse Jacobian Function Jac_SP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 USE CB05TUCl_EPA_Model

      INTEGER nlhs, nrhs
      INTEGER plhs(*), prhs(*)
      INTEGER mxGetPr, mxCreateFull, mxGetM, mxgetN
      INTEGER VPtr, FPtr, RPtr, JVSPtr
      REAL(kind=dp) V(82), F(4), RCT(205)
      REAL(kind=dp) JVS(861)

! Check for the right number of input arguments
      IF ( nrhs .ne. 3 ) THEN
         CALL mexErrMsgTxt('Jac_SP requires 3 input vectors: &
     &V(82), F(4), RCT(205)')
      END IF 
! Check for the right number of output arguments
      IF ( nlhs .ne. 1 ) THEN
         CALL mexErrMsgTxt('Jac_SP requires 1 output vector: &
     &JVS(861)')
      END IF 

      plhs(1) = mxCreateDoubleMatrix(861,1,0)

      VPtr = mxGetPr(prhs(1))
      CALL mxCopyPtrToReal8(VPtr,V,82)
      
      FPtr = mxGetPr(prhs(2))
      CALL mxCopyPtrToReal8(FPtr,F,4)
      
      RPtr = mxGetPr(prhs(3))
      CALL mxCopyPtrToReal8(RPtr,RCT,205)

      JVSPtr = mxGetPr(plhs(1))

      CALL Jac_SP( V, F, RCT, JVS )

      CALL mxCopyReal8ToPtr(JVS, JVSPtr, 861)

 END SUBROUTINE mexFunction

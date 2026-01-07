!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Matlab Gateway for the Sparse Jacobian Function Jac_SP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 USE MOZART_T1_Model

      INTEGER nlhs, nrhs
      INTEGER plhs(*), prhs(*)
      INTEGER mxGetPr, mxCreateFull, mxGetM, mxgetN
      INTEGER VPtr, FPtr, RPtr, JVSPtr
      REAL(kind=dp) V(155), F(5), RCT(360)
      REAL(kind=dp) JVS(1826)

! Check for the right number of input arguments
      IF ( nrhs .ne. 3 ) THEN
         CALL mexErrMsgTxt('Jac_SP requires 3 input vectors: &
     &V(155), F(5), RCT(360)')
      END IF 
! Check for the right number of output arguments
      IF ( nlhs .ne. 1 ) THEN
         CALL mexErrMsgTxt('Jac_SP requires 1 output vector: &
     &JVS(1826)')
      END IF 

      plhs(1) = mxCreateDoubleMatrix(1826,1,0)

      VPtr = mxGetPr(prhs(1))
      CALL mxCopyPtrToReal8(VPtr,V,155)
      
      FPtr = mxGetPr(prhs(2))
      CALL mxCopyPtrToReal8(FPtr,F,5)
      
      RPtr = mxGetPr(prhs(3))
      CALL mxCopyPtrToReal8(RPtr,RCT,360)

      JVSPtr = mxGetPr(plhs(1))

      CALL Jac_SP( V, F, RCT, JVS )

      CALL mxCopyReal8ToPtr(JVS, JVSPtr, 1826)

 END SUBROUTINE mexFunction

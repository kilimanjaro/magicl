;;;; tests/matrix-tests.lisp
;;;;
;;;; Author: Cole Scott

(in-package #:magicl-tests)

(defmacro is-matrix (&rest rest)
  `(progn ,@(loop :for m :in rest
                  :collect `(is (subtypep (type-of ,m) 'matrix)))))

(defmacro is-not-matrix (&rest rest)
  `(progn ,@(loop :for m :in rest
                  :collect `(is (not (subtypep (type-of ,m) 'matrix))))))

(deftest test-identity-matrix-p ()
  "Test that identity matrices can be identified by IDENTITY-MATRIX-P for all types of matrixes from 1x1 to 64x64"
  (dolist (type +magicl-types+)
    (loop :for i :from 1 :to 64 :do
      (is (magicl:identity-matrix-p (magicl:eye (list i i) :type type)))
      (is (not (magicl:identity-matrix-p (magicl:eye (list i i) :value 2 :type type))))
      (is (not (magicl:identity-matrix-p (magicl:eye (list i i) :offset 1 :type type))))
      (is (not (magicl:identity-matrix-p (magicl:const 0 (list i i) :type type)))))))

(deftest test-square-matrix-p ()
  "Test that square matrices can be identified by IDENTITY-MATRIX-P for all types of matrixes from 1x1 to 64x64"
  (dolist (type +magicl-types+)
    (loop :for i :from 1 :to 64 :do
      (is (magicl:square-matrix-p (magicl:empty (list i i) :type type)))
      (is (not (magicl:square-matrix-p (magicl:empty (list i (* 2 i)) :type type)))))))

;; Multiplication

(deftest test-matrix-multiplication ()
  "Test multiplication for random pairs of matrices"
  (labels ((mult (a b)
             (let* ((m (magicl:nrows a))
                    (n (magicl:ncols b))
                    (p (magicl:ncols a))
                    (target (magicl:empty (list m n))))
               (loop :for i :below m :do
                 (loop :for j :below n :do
                   (setf (magicl:tref target i j)
                         (loop :for k :below p
                               :sum (* (magicl:tref a i k)
                                       (magicl:tref b k j))))))
               target)))
    (dolist (magicl::*default-tensor-type* +magicl-float-types+)
      (loop :for i :below 1000 :do
        (let* ((n (1+ (random 5)))
               (m (1+ (random 5)))
               (k (1+ (random 5)))
               (a (magicl:rand (list m k)))
               (b (magicl:rand (list k n)))
               (c (magicl:rand (list m n))))
          ;; Check that multiplication returns the correct result
          (is (magicl:=
               (mult a b)
               (magicl:mult a b)))
          ;; Check that transposing doesn't affect correctness
          (is (magicl:=
               (mult (magicl:transpose a) c)
               (magicl:mult a c :transa :t)))
          (is (magicl:=
               (mult b (magicl:transpose c))
               (magicl:mult b c :transb :t)))
          (is (magicl:=
               (mult (magicl:transpose b) (magicl:transpose a))
               (magicl:mult b a :transa :t :transb :t)))
          ;; Check that alpha correctly scales the matrices
          (is (magicl:=
               (mult (magicl:scale a 2) b)
               (magicl:mult a b :alpha (coerce 2 magicl::*default-tensor-type*)))))))))

(deftest test-matrix-vector-multiplication ()
  "Test multiplication for random pairs of matrix and vectors"
  (labels ((mult (a x)
             (assert (= (magicl:ncols a) (magicl:size x)))
             (let* ((m (magicl:nrows a))
                    (n (magicl:ncols a))
                    (target (magicl:empty (list m))))
               (loop :for i :below m :do
                 (setf (magicl:tref target i)
                       (loop :for k :below n
                             :sum (* (magicl:tref a i k)
                                     (magicl:tref x k)))))
               target)))
    (dolist (magicl::*default-tensor-type* +magicl-float-types+)
      (loop :for i :below 1000 :do
        (let* ((n (1+ (random 5)))
               (m (1+ (random 5)))
               (a (magicl:rand (list m n)))
               (x (magicl:rand (list n)))
               (y (magicl:rand (list m))))

          ;; Check that multiplication returns the correct result
          (is (magicl:=
               (mult a x)
               (magicl:mult a x)))

          ;; Check that transposing doesn't affect correctness
          (is (magicl:=
               (mult (magicl:transpose a) y)
               (magicl:mult a y :transa :t)))

          ;; Check that alpha correctly scales the matrices
          (is (magicl:=
               (mult (magicl:scale a 2) x)
               (magicl:mult a x :alpha (coerce 2 magicl::*default-tensor-type*)))))))))

(deftest test-matrix-multiplication-errors ()
  (signals simple-error (magicl:@
                         (magicl:empty '(3 3))
                         (magicl:empty '(1 1))))
  (signals simple-error (magicl:@
                         (magicl:empty '(1 2))
                         (magicl:empty '(1 2))))
  (signals simple-error (magicl:@
                         (magicl:empty '(5 2))
                         (magicl:empty '(2 3))
                         (magicl:empty '(2 3)))))

(deftest test-complex-matrix-multiplication-results ()
  "Test a few basic complex matrix multiplications"
  (let* ((m (magicl:from-list '(#C(1d0 2d0) #C(3d0 4d0) #C(5d0 6d0) #C(7d0 8d0)) '(2 2) :layout :row-major))
         (m-old (magicl::deep-copy-tensor m))
         (x (magicl:from-list '(#C(1d0 2d0) #C(3d0 4d0)) '(2 1)))
         (x-old (magicl::deep-copy-tensor x))
         (expected (magicl:from-list '(#C(-10d0 28d0) #C(-18d0 68d0)) '(2 1))))
    ;; Check that the multiplication is correct and does not throw any errors
    (is (magicl:= expected (magicl:@ m x)))

    ;; Check that the multiplication did not modify the inputs
    (is (magicl:= m-old m))
    (is (magicl:= x-old x))

    ;; Check that doing 2x1 @ 2x2 errors
    (signals error (magicl:@ x m))))

;;; Block Matrix Routines

(deftest test-block-diagonal ()
  "Test that we can construct block diagonal matrices."
  (let ((expected (magicl:from-list '(0d0 0d0 0d0
                                      0d0 1d0 1d0
                                      0d0 1d0 1d0)
                                    '(3 3))))
    (is (magicl:= expected
                  (magicl:block-diag
                   (list (magicl:zeros '(1 1))
                         (magicl:ones '(2 2))))))))

(deftest test-matrix-stacking ()
  "Test that we can stack matrices 'horizontally' and 'vertically'."
  (let ((expected (magicl:from-list '(1 2 3
                                      4 5 6)
                                    '(2 3))))
    (is (magicl:= expected
                  (magicl:hstack
                   (loop :for j :below 3
                         :collect (magicl:column expected j)))))
    (is (magicl:= expected
                  (magicl:vstack
                   (loop :for i :below 2
                         :collect (magicl:row expected i)))))))


(deftest test-block-matrix-construction ()
  "Test that we can construct a block matrix."
  (let ((mat
          (magicl::block-matrix (list (magicl:zeros '(3 2))     (magicl:eye 3 :value 3d0)
                                      (magicl:eye 2 :value 2d0)     (magicl:zeros '(2 3)))
                                '(2 2))))
    (is (magicl:=
         mat
         (magicl:from-list '(0d0 0d0 3d0 0d0 0d0
                             0d0 0d0 0d0 3d0 0d0
                             0d0 0d0 0d0 0d0 3d0
                             2d0 0d0 0d0 0d0 0d0
                             0d0 2d0 0d0 0d0 0d0)
                           '(5 5))))))


;; QR and Eigenvalue Tests

(deftest test-qr-special-cases ()
  "Test that the QR factorization works as advertised in a few silly cases."
  (dolist (mat (list
		(magicl:eye 3)
		(magicl:ones '(3 3))
		(magicl:ones '(3 1))
                (magicl:ones '(1 3))
		(magicl:from-list '(#C(1d0 0d0) #C(0d0 1d0) #C(1d0 1d0) #C(0d0 0d0)) '(2 2))))
    (multiple-value-bind (Q R) (magicl:qr mat)
      (is (magicl:= (magicl:@ (magicl:dagger Q) Q)
		    (magicl:eye (magicl:ncols Q) :type (magicl:element-type Q))))
      (is (magicl:= R (magicl:upper-triangular R)))
      (is (magicl:= mat (magicl:@ Q R))))))

(deftest test-hermitian-eig ()
  "Test that we can compute eigenvectors & values of Hermitian matrices."
  (let ((matrix-size 5)
        (repetitions 5)
        (eig-range 1d0)
        (element-types (list 'double-float '(complex double-float))))
    (dolist (element-type element-types)
      (dotimes (i repetitions)
        (let* ((evals (sort (loop :for i :below matrix-size
                                  :collect (random eig-range))
                            #'<))
               (D (magicl:from-diag evals :type element-type))
               (Q (magicl:random-unitary (list matrix-size matrix-size) :type element-type))
               (M (magicl:@ Q D (magicl:conjugate-transpose Q))))
          (is (magicl:identity-matrix-p (magicl:@ Q (magicl:conjugate-transpose Q))))
          (is (magicl:hermitian-matrix-p m))
          (multiple-value-bind (%evals Q)
              (magicl::hermitian-eig m)
            (is (every (lambda (a b) (< (abs (- a b)) 1d-8))
                       evals
                       (sort (mapcar #'realpart %evals) #'<)))
            (is (magicl:= m (magicl:@ Q (magicl:from-diag %evals :type element-type) (magicl:dagger Q))))))))))

(deftest test-bidiagonal ()
  "Test that we reduce an arbitrary matrix to bidiagonal form."
  (let* ((u (magicl:random-unitary (list 5 5)))
         (d (magicl:from-diag '(1d0 2d0 3d0 4d0 5d0)))
         (v (magicl:random-unitary (list 5 5)))
         (m (magicl:@ u d (magicl:dagger v))))
    (multiple-value-bind (u b v) (magicl::bidiagonal m)
      (is (magicl:= m (magicl:@ u b (magicl:dagger v)))))))

(deftest test-svd ()
  "Test the full and reduced SVDs."
  (labels ((mul-diag-times-gen (diag matrix)
             "Returns a newly allocated matrix resulting from the product of DIAG (a diagonal real matrix) with MATRIX (a complex matrix)."
             #+ignore
             (declare (type matrix diag matrix)
                      (values matrix))
             (let* ((m (magicl:nrows diag))
                    (k (magicl:ncols matrix))
                    (result (magicl:empty (list m k))))
               (dotimes (i (min m (magicl:ncols diag)) result)
                 (let ((dii (magicl:tref diag i i)))
                   (dotimes (j k)
                     (setf (magicl:tref result i j)
                           (* dii (magicl:tref matrix i j))))))))

           (norm-inf (matrix)
             "Return the infinity norm of vec(MATRIX)."
             (let ((data (magicl::storage matrix)))
               (reduce #'max data :key #'abs)))

           (zero-p (matrix &optional (tolerance (* 1.0d2 double-float-epsilon)))
             "Return T if MATRIX is close to zero (within TOLERANCE)."
             (< (norm-inf matrix) tolerance))

           (check-full-svd (matrix)
             "Validate full SVD of MATRIX."
             (let ((m (magicl:nrows matrix))
                   (n (magicl:ncols matrix)))
               (multiple-value-bind (u sigma vh)
                   (magicl:svd matrix)
                 (is (= (magicl:nrows u) (magicl:ncols u) m))
                 (is (and (= (magicl:nrows sigma) m) (= (magicl:ncols sigma) n)))
                 (is (= (magicl:nrows vh) (magicl:ncols vh) n))
                 (print 

(deftest test-hermitian-eig ()
  "Test that we can compute eigenvectors & values of Hermitian matrices."
  (let ((matrix-size 5)
        (repetitions 5)
        (eig-range 1d0)
        (element-types (list 'double-float '(complex double-float))))
    (dolist (element-type element-types)
      (dotimes (i repetitions)
        (let* ((evals (sort (loop :for i :below matrix-size
                                  :collect (random eig-range))
                            #'<))
               (D (magicl:from-diag evals :type element-type))
               (Q (magicl:random-unitary (list matrix-size matrix-size) :type element-type))
               (M (magicl:@ Q D (magicl:conjugate-transpose Q))))
          (is (magicl:identity-matrix-p (magicl:@ Q (magicl:conjugate-transpose Q))))
          (is (magicl:hermitian-matrix-p m))
          (multiple-value-bind (%evals Q)
              (magicl::hermitian-eig m)
            (is (every (lambda (a b) (< (abs (- a b)) 1d-8))
                       evals
                       (sort (mapcar #'realpart %evals) #'<)))
            (is (magicl:= m (magicl:@ Q (magicl:from-diag %evals :type element-type) (magicl:dagger Q))))))))))

(deftest test-bidiagonal ()
  "Test that we reduce an arbitrary matrix to bidiagonal form."
  (let* ((u (magicl:random-unitary (list 5 5)))
         (d (magicl:from-diag '(1d0 2d0 3d0 4d0 5d0)))
         (v (magicl:random-unitary (list 5 5)))
         (m (magicl:@ u d (magicl:dagger v))))
    (multiple-value-bind (u b v) (magicl::bidiagonal m)
      (is (magicl:= m (magicl:@ u b (magicl:dagger v)))))))

(deftest test-svd ()
  "Test the full and reduced SVDs."
  (labels ((mul-diag-times-gen (diag matrix)
             "Returns a newly allocated matrix resulting from the product of DIAG (a diagonal real matrix) with MATRIX (a complex matrix)."
             #+ignore
             (declare (type matrix diag matrix)
                      (values matrix))
             (let* ((m (magicl:nrows diag))
                    (k (magicl:ncols matrix))
                    (result (magicl:empty (list m k))))
               (dotimes (i (min m (magicl:ncols diag)) result)
                 (let ((dii (magicl:tref diag i i)))
                   (dotimes (j k)
                     (setf (magicl:tref result i j)
                           (* dii (magicl:tref matrix i j))))))))

           (norm-inf (matrix)
             "Return the infinity norm of vec(MATRIX)."
             (let ((data (magicl::storage matrix)))
               (reduce #'max data :key #'abs)))

           (zero-p (matrix &optional (tolerance (* 1.0d2 double-float-epsilon)))
             "Return T if MATRIX is close to zero (within TOLERANCE)."
             (< (norm-inf matrix) tolerance))

           (check-full-svd (matrix)
             "Validate full SVD of MATRIX."
             (let ((m (magicl:nrows matrix))
                   (n (magicl:ncols matrix)))
               (multiple-value-bind (u sigma vh)
                   (magicl:svd matrix)
                 (is (= (magicl:nrows u) (magicl:ncols u) m))
                 (is (and (= (magicl:nrows sigma) m) (= (magicl:ncols sigma) n)))
                 (is (= (magicl:nrows vh) (magicl:ncols vh) n))
                 (is (zero-p (magicl:.- matrix (magicl:@ u (mul-diag-times-gen sigma vh))))))))

           (check-reduced-svd (matrix)
             "Validate reduced SVD of MATRIX."
             (let* ((m (magicl:nrows matrix))
                    (n (magicl:ncols matrix))
                    (k (min m n)))

               (multiple-value-bind (u sigma vh)
                   (magicl:svd matrix :reduced t)
                 (is (and (= (magicl:nrows u) m)
                          (= (magicl:ncols u) k)))
                 (is (= (magicl:nrows sigma) (magicl:ncols sigma) k))
                 (is (and (= (magicl:nrows vh) k)
                          (= (magicl:ncols vh) n)))
                 (is (zero-p (magicl:.- matrix (magicl:@ u (mul-diag-times-gen sigma vh)))))))))

    (let ((tall-thin-matrix (magicl:rand '(8 2))))
      (check-full-svd tall-thin-matrix)
      (check-reduced-svd tall-thin-matrix))

    (let ((short-fat-matrix (magicl:rand '(2 8))))
;      (check-full-svd short-fat-matrix)
;     (check-reduced-svd short-fat-matrix)
      )))

;;;; double-float.lisp
;;;;
;;;; Author: Cole Scott

(in-package #:magicl)

(deftensor tensor/int32 (signed-byte 32))

(defmatrix matrix/int32 (signed-byte 32) tensor/int32)

(defvector vector/int32 (signed-byte 32) tensor/int32)

(defcompatible
    (lambda (tensor)
      (case (order tensor)
        (1 '(vector/int32
             tensor/int32))
        (2 '(matrix/int32
             tensor/int32))
        (t '(tensor/int32))))
  tensor/int32
  matrix/int32
  vector/int32)

(register-matrix-matrix-multiply matrix/int32 (signed-byte 32))
(register-matrix-vector-multiply matrix/int32 vector/int32 (signed-byte 32))

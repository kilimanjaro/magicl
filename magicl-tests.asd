;;;; magicl-tests.asd
;;;;
;;;; Author: Joseph Lin
;;;;         Erik Davis

(asdf:defsystem #:magicl-tests/core
  :license "BSD 3-Clause (See LICENSE.txt)"
  :description "Regression tests for MAGICL."
  :author "Rigetti Computing"
  :depends-on (#:alexandria
               #:uiop
               #:magicl/core
               #:fiasco)
  :perform (asdf:test-op (o s)
                         (uiop:symbol-call :magicl-tests
                                           '#:run-magicl-tests))
  :pathname "tests/"
  :serial t
  :components ((:file "packages")
               (:file "suite")
               (:file "constants")
               (:file "util-tests")
               (:file "abstract-tensor-tests")
               (:file "specialization-tests")
               (:file "constructor-tests")
               (:file "matrix-tests")))

(asdf:defsystem #:magicl-tests
  :license "BSD 3-Clause (See LICENSE.txt)"
  :description "Regression tests for MAGICL."
  :author "Rigetti Computing"
  :depends-on (#:alexandria
               #:uiop
               #:magicl
               #:magicl-examples
               #:fiasco)
  :perform (asdf:test-op (o s)
                         (uiop:symbol-call :magicl-tests
                                           '#:run-magicl-tests))
  :pathname "tests/"
  :serial t
  :components ((:file "packages")
               (:file "suite")
               (:file "constants")
               (:file "util-tests")
               (:file "abstract-tensor-tests")
               (:file "specialization-tests")
               (:file "constructor-tests")
               (:file "matrix-tests")
               (:module "extensions"
                :serial t
                :components ((:file "high-level-tests")))))

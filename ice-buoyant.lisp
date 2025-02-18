;; (restrict-compiler-policy 'speed 3 3)
;; (restrict-compiler-policy 'debug 0 0)
;; (restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;; (ql:quickload "cl-mpm/examples/ice-visco")
(in-package :cl-mpm/examples/ice-buoyancy)

(defun setup (&key (refine 1) (mps 2))
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let* ((density 916.7d0)
           (mesh-resolution (/ 10d0 refine))
           (offset (* mesh-resolution 2))
           (datum (+ 100d0 offset))
           (ice-height 200)
           (aspect 1)
           (domain-size (list 800d0 400d0))
           (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
           (block-size (list (* ice-height aspect) ice-height)))
      (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                                 :sim-type
                                                 'cl-mpm/damage::mpm-sim-damage
                                                 ;; 'cl-mpm::mpm-sim-usf
                                                 ))
      (let* ((E 1d9)
             (angle 40d0)
             (init-c 1d5)
             (init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180))))
             (gf 1000d0)
             (length-scale (* mesh-resolution 2d0))
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
             (oversize (cl-mpm/damage::compute-oversize-factor (- 1d0 1d-6) ductility)))
        (defparameter *removal-strain* (* 100d0 (/ init-stress 1)))
        (when (= rank 0)
          (format t "Mesh size ~F~%" mesh-resolution)
          (format t "Estimated oversize ~F~%" oversize)
          (format t "Estimated lc ~E~%" length-scale)
          (format t "Estimated ductility ~E~%" ductility)
          (format t "Init stress ~E~%" init-stress)
          (format t "Removal strain ~E~%" *removal-strain*))
        (cl-mpm:add-mps
         *sim*
         (cl-mpm/setup:make-block-mps
          (list 0 offset)
          block-size
          (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
          density
          'cl-mpm/particle::particle-chalk-delayed
          :E 1d9
          :nu 0.24d0

          :kt-res-ratio 1d0
          :kc-res-ratio 0d0
          :g-res-ratio 0.1d0

          :initiation-stress init-stress;18d3
          :friction-angle angle
          :psi (* 0d0 (/ pi 180))
          :phi (* angle (/ pi 180))
          :c (* init-c oversize)
          :softening 0d0
          :ductility ductility
          :local-length length-scale
          :delay-time 1d1
          :delay-exponent 2
          :enable-plasticity nil
          :enable-damage t
          ;; 'cl-mpm/particle::particle-finite-viscoelastic-ice
          ;; :E 1d9
          ;; :nu 0.325d0
          ;; :visc-factor 11.1d6
          ;; :visc-power 3d0
          ;; :enable-viscosity nil

          :gravity -9.8d0
          ))
        ;; (cl-mpm:add-mps
        ;;  *sim*
        ;;  (cl-mpm/setup:make-block-mps
        ;;   (list (first block-size) offset)
        ;;   block-size
        ;;   (mapcar (lambda (e) (* (/ e mesh-resolution) mps)) block-size)
        ;;   density
        ;;   'cl-mpm/particle::particle-fluid
        ;;   :rest-density 1000d0
        ;;   :stiffness 1d3
        ;;   :adiabatic-index 7
        ;;   :viscosity 1d-3
        ;;   ;; :E 1d9
        ;;   ;; :nu 0.24d0
        ;;   :gravity -9.8d0
        ;;   ))
        )
      (setf
       (cl-mpm:sim-bcs *sim*)
       (cl-mpm/bc::make-outside-bc-varfix
        (cl-mpm:sim-mesh *sim*)
        '(0 nil 0)
        '(0 nil 0)
        '(nil 0 0)
        '(nil 0 0)
        '(nil nil 0)
        '(nil nil 0)))

      (setf (cl-mpm:sim-mass-scale *sim*) 1d4)
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.1d0
               (sqrt 1d4)
               (cl-mpm/setup:estimate-critical-damping *sim*)))
      ;; (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-4)
      (setf (cl-mpm::sim-enable-fbar *sim*) t)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (setf (cl-mpm::sim-allow-mp-split *sim*) t)
      ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
      (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)



      (setf (cl-mpm:sim-dt *sim*)
            (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
      (setf (cl-mpm::sim-enable-damage *sim*) nil)
      (setf *run-sim* t)
      (if t
          (cl-mpm:add-bcs-force-list
           *sim*
           (cl-mpm/buoyancy::make-bc-buoyancy-clip
            *sim*
            datum
            1000d0
            (lambda (pos datum) t)
            :visc-damping (* 1d-3 (cl-mpm/setup:estimate-critical-damping *sim*))
            ))
          (cl-mpm:add-bcs-force-list
           *sim*
           (cl-mpm/buoyancy::make-bc-buoyancy-body
            *sim*
            datum
            1000d0
            (lambda (pos) t))))
      (let ((domain-half (* 0.5d0 (first domain-size)))
            (friction 0d0))
        (defparameter *floor-bc*
          (cl-mpm/penalty::make-bc-penalty-distance-point
           *sim*
           (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           (cl-mpm/utils:vector-from-list (list
                                           domain-half
                                           offset
                                           0d0))
           (* domain-half 1.1d0)
           (* 1d9 0.1d0)
           friction
           0d0)))

      ;; (defparameter *bc-erode*
      ;;   (cl-mpm/erosion::make-bc-erode
      ;;    *sim*
      ;;    :enable nil
      ;;    :rate 1d-1
      ;;    :scalar-func (lambda (pos)
      ;;                   (min 1d0 (exp (* 0.5d0 (- (cl-mpm/utils:varef pos 1) datum)))))
      ;;    :clip-func (lambda (pos)
      ;;                 (>= datum (cl-mpm/utils:varef pos 1))
      ;;                 ;; (and (>= (+ offset mesh-resolution) (cl-mpm/utils:varef pos 1)))
      ;;                 )))
      (cl-mpm:add-bcs-force-list
       *sim*
       *floor-bc*
       )
      ;; (cl-mpm:add-bcs-force-list
      ;;  *sim*
      ;;  *bc-erode*
      ;;  )
      (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
      )))

(defun run (&key (output-dir "./output/"))
  (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
  ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let ((dt-scale 0.5d0))
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :oobf-crit 1d-1
       :energy-crit 1d-1
       :conv-steps 1000
       :dt-scale dt-scale
       :post-iter-step
       (lambda (i oobf energy)
         (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_conv_~2,'0d_~5,'0d.vtk" rank i)) *sim*)
         )))

    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (when (typep mp 'cl-mpm/particle::particle-finite-viscoelastic-ice)
         (setf (cl-mpm/particle::mp-enable-viscosity mp) t))))

    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d8)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-4
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup:estimate-critical-damping *sim*)))

    ;; (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
    (let* ((dt-scale 0.5d0)
           (target-time 1d4)
           (work 0d0)
           (oobf 0d0)
           (energy 0d0))
      (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
      (setf substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
      (when (= rank 0)
        (format t "Substeps ~D~%" substeps))

      (loop for step from 0 below 1000
            while *run-sim*
            do
            (progn
              (let ((mp-arr (make-array (cl-mpi:mpi-comm-size) :element-type 'double-float :initial-element 0d0)))
                (setf (aref mp-arr rank) (float (length (cl-mpm:sim-mps *sim*)) 0d0))
                (cl-mpm/mpi::mpi-vector-sum mp-arr)
                (when (= rank 0)
                  (format t "Step ~D~%" step)
                  (format t "Step ~D~%" step)
                  (format t "MPs: ~{~D ~}" (mapcar #'round (map 'list #'identity mp-arr)))
                  ))
              (cl-mpm/mpi::load-balance-algo *sim* 
                                     :step-size 1d-2
                                     ;:min-bounds 1.01d0
                                     ;:max-bounds 1.2d0
                                     ;:max-bounds (* 2d0 *balance-point* )
                                     )
              (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_~2,'0d_~5,'0d.vtk" rank step)) *sim*)
              (cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank step)) *sim*)
              (setf work 0d0)
              (time
                (dotimes (i substeps)
                  (cl-mpm:update-sim *sim*)
                  ;; (cl-mpm::remove-mps-func
                  ;;  *sim*
                  ;;  (lambda (mp)
                  ;;    (> (magicl:det (cl-mpm/particle::mp-deformation-gradient mp)) 1.2d0)))
                  ;; (cl-mpm::split-mps-eigenvalue *sim*)
                  (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))))
              (setf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
              (setf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
              (setf
                energy (/ energy substeps)
                oobf (/ oobf substeps))
              (if (= work 0d0)
                  (setf energy 0d0)
                  (setf energy (abs (/ energy work))))
              (when (= rank 0)
                (format t "OOBF ~E - Energy ~E~%" oobf energy)))
            ))))

(defun run-adaptive (&key (output-dir "./output/")
                       )
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
    ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    (defparameter *damping-water* 0d0)
    (let ((dt-scale 0.5d0))
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :oobf-crit 1d-1
       :energy-crit 1d-1
       :dt-scale dt-scale
       :post-iter-step
       (lambda (i oobf energy)
         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~2,'0d_~5,'0d.vtk" rank i)) *sim*)
         )))

    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 0d-6
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    ;; (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
    (let* ((dt-scale 0.5d0)
           (substeps 0d0)
           (work 0d0)
           (oobf 0d0)
           (energy 0d0)
           (sim-state :accelerate)
           (accelerate-target-time 1d1)
           (accelerate-mass-scale 1d4)
           (collapse-target-time 0.1d0)
           (collapse-mass-scale 1d0)
           (criteria-energy 1d-2)
           (criteria-oobf 1d-1)
           (target-time 1d0)
           (time 0d0)
           )
      (setf (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale
            target-time accelerate-target-time)
      (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
      (setf substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
      (format t "Substeps ~D~%" substeps)
      (loop for step from 0 below 1000
            while *run-sim*
            do
               (when (= rank 0)
                 (format t "Step ~D~%" step))
               ;; (setf work 0d0)
               (cl-mpm/mpi::load-balance-algo *sim* :step-size 1d-2)
               (setf energy 0d0
                     oobf 0d0
                     )
               (time
                (dotimes (i substeps)
                  (cl-mpm:update-sim *sim*)
                  (incf time (cl-mpm:sim-dt *sim*))

                  ;; (cl-mpm::remove-mps-func
                  ;;  *sim*
                  ;;  (lambda (mp)
                  ;;    ;; (> (cl-mpm/fastmaths:det (cl-mpm/particle:mp-deformation-gradient mp)) 1.5d0)
                  ;;    (multiple-value-bind (s1 s2 s3) (cl-mpm/utils::principal-stresses-3d (cl-mpm/particle::mp-undamaged-stress mp) )
                  ;;      (> (max s1 s2 s3) *removal-strain*))
                  ;;    ))
                  ;; (cl-mpm::split-mps-eigenvalue *sim*)
                  (incf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                  (incf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                  (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))))
               (setf
                energy (/ energy substeps)
                oobf (/ oobf substeps))

               (if (= work 0d0)
                   (setf energy 0d0)
                   (setf energy (abs (/ energy work))))

               (if (or
                    (> energy criteria-energy)
                    (> oobf criteria-oobf)
                    ;; t
                    ;; nil
                    ;; (> work 1d6)
                    )
                   (when (not (eq sim-state :collapse))
                     (setf sim-state :collapse)
                     (when (= rank 0)
                       (format t "Changed to collapse~%"))
                     (setf work 0d0)
                     )
                   (progn
                     (when (not (eq sim-state :accelerate))
                       (when (= rank 0)
                         (format t "Changed to accelerate~%"))
                       (setf work 0d0)
                       (setf sim-state :accelerate)
                       ;; (cl-mpm::remove-mps-func
                       ;;  *sim*
                       ;;  (lambda (p)
                       ;;    (and (> (cl-mpm::mp-damage p) 0.99d0)
                       ;;         (= (cl-mpm/particle::mp-index p) 0))
                       ;;    ))
                       (cl-mpm:iterate-over-mps
                        (cl-mpm:sim-mps *sim*)
                        (lambda (mp)
                          ;; (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-acceleration mp))
                          (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))
                          ))
                       )))
               (case sim-state
                 (:accelerate
                  (when (= rank 0)
                    (format t "Accelerate timestep~%"))
                  (setf
                   target-time accelerate-target-time
                   (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale))
                 (:collapse
                  (when (= rank 0)
                    (format t "Collapse timestep~%"))
                  (setf
                   ;; work 0d0
                   target-time collapse-target-time
                   (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale)))


               (when (= rank 0)
                 (format t "OOBF ~E - Energy ~E~%" oobf energy)
                 (format t "State ~A~%" sim-state))

               (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~2,'0d_~5,'0d.vtk" rank step)) *sim*)
               (cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank step)) *sim*)

               (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                 (when (= rank 0)
                   (format t "CFL dt estimate: ~f~%" dt-e)
                   (format t "CFL step count estimate: ~D~%" substeps-e))
                 (setf substeps substeps-e))
               (swank.live:update-swank)))))

(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (setup :refine 2 :mps 2)
    (when (typep *sim* 'cl-mpm/mpi::mpm-sim-mpi)
   	  (let* ((height 1)
             (dsize (ceiling (cl-mpi:mpi-comm-size) height)))
        (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize height 1)))
      (cl-mpm/mpi::setup-domain-bounds *sim*)
      (setf cl-mpm/mpi::*prune-nodes* nil)
      (cl-mpm/mpi::load-balance-algo *sim* 
                                     :step-size 1d-2
                                        ;:min-bounds 1.01d0
                                        ;:max-bounds 1.2d0
                                        ;:max-bounds (* 2d0 *balance-point* )
                                     )
      (cl-mpm/mpi::domain-decompose *sim*))

    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps *sim*)))
    (when (= rank 0)
      ;; (pprint (mpm-sim-mpi-domain-bounds *sim*))
      (format t "Run mpi~%"))
    (run-adaptive :output-dir *output-directory*)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  )

;; (defparameter *output-directory* (merge-pathnames "/nobackup/rmvn14/paper-2/ice-buoyant/"))
(defparameter *output-directory* (merge-pathnames "./output/"))
(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "16"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
;; (setup)
;; (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;; (run :output-dir *output-directory*)



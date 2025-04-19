;; (restrict-compiler-policy 'speed 3 3)
;; (restrict-compiler-policy 'debug 0 0)
;; (restrict-compiler-policy 'safety 0 0)
;; (setf *block-compile-default* t)
;; (ql:quickload "cl-mpm/examples/ice-visco")
(in-package :cl-mpm/examples/ice-buoyancy)

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage) 
  (declare (double-float local-length damage)) 
  (* local-length (max (sqrt (- 1d0 damage)) 1d-10)))

(defmethod cl-mpm::update-particle (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt)
  (cl-mpm::update-particle-kirchoff mesh mp dt)
  ;(cl-mpm::update-domain-det mesh mp dt)
  ;(cl-mpm::update-domain-polar-2d mesh mp dt)
  ;(cl-mpm::co-domain-corner-2d mesh mp dt)
  (cl-mpm::update-domain-polar-2d mesh mp dt)
  ;(cl-mpm::update-domain-midpoint mesh mp dt)
  ;(cl-mpm::update-domain-deformation mesh mp dt)
  ;(cl-mpm::scale-domain-size mesh mp)
  )

(defun setup (&key (refine 1) 
                   (mps 2)
                   (height 800)
                   (cutback 0d0)
                   (aspect 1)
                   (fbar nil)
                )
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let* ((density 918d0)
           (water-density 1028)
           (mesh-resolution (/ 10d0 refine))
           (offset (* mesh-resolution 0))
           (ice-height height)
           (floating-point (* ice-height (/ density water-density)))
           (water-level (* floating-point 0.8d0))
           (datum (* (round (+ water-level offset) mesh-resolution) mesh-resolution))
           (ice-length (float (* aspect ice-height) 0d0))
           (domain-size (list (* 2 ice-length) (* 2 ice-height)))
           (element-count (mapcar (lambda (x) (round x mesh-resolution)) domain-size))
           (block-size (list (* ice-height aspect) ice-height))
           )
      (defparameter *water-height* datum)
      (defparameter *ice-height* (+ ice-height offset))
      (defparameter *ice-length* ice-length)
      (setf *sim* (cl-mpm/setup::make-simple-sim mesh-resolution element-count
                                                 :sim-type
                                                 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
                                                 ;'cl-mpm/damage::mpm-sim-damage
                                                 ;; 'cl-mpm::mpm-sim-usf
                                                 ))
      (let* ((E 1d9)
             ;(angle 40d0)
             ;(init-c 1d5)
             ;(init-stress (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile init-c (* angle (/ pi 180)))) 
             (angle 40d0) 
             (init-stress (* 0.1185d6 1d0)) 
             (init-c (cl-mpm/damage::mohr-coloumb-tensile-to-coheasion init-stress (* angle (/ pi 180))))
             (gf 1000d0)
             ;(length-scale (* mesh-resolution 2d0))
             (length-scale 20d0)
             (delay-time 1d3)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress E))
             (oversize (cl-mpm/damage::compute-oversize-factor (- 1d0 1d-3) ductility)))
        (defparameter *delay-time* delay-time)
        (defparameter *removal-strain* (* 100d0 (/ init-stress 1)))
        (when (= rank 0)
          (format t "Mesh size ~F~%" mesh-resolution)
          (format t "Ice height ~F~%" ice-height)
          (format t "Cliff height ~F~%" (- ice-height datum))
          (format t "Estimated oversize ~F~%" oversize)
          (format t "Estimated lc ~E~%" length-scale)
          (format t "Estimated ductility ~E~%" ductility)
          (format t "Init stress ~E~%" init-stress)
          (format t "Init c ~E~%" init-c)
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
          :g-res-ratio 0.5d0

          :initiation-stress init-stress;18d3
          :friction-angle angle
          :psi (* 0d0 (/ pi 180))
          :phi (* angle (/ pi 180))
          :c (* init-c oversize)
          :softening 0d0
          :ductility ductility
          :local-length length-scale
          :delay-time delay-time
          :delay-exponent 2
          :enable-plasticity nil
          :enable-damage t
          :gravity -9.8d0
          )))
      ;; (let ((k 1d0))
      ;;   (cl-mpm/setup::initialise-stress-self-weight-vardatum
      ;;    *sim*
      ;;    (lambda (pos)
      ;;      (let ((alpha (/ (- (cl-mpm/utils::varef pos 0)
      ;;                         ice-length) ice-length)))
      ;;        (+ offset
      ;;           (* alpha end-height)
      ;;           (* (- 1d0 alpha) start-height))))
      ;;    ;; :k-x k
      ;;    ;; :k-z k
      ;;    ))
      (cl-mpm/setup::initialise-stress-self-weight *sim* (+ offset ice-height))
      (when (> cutback 0d0)
        (cl-mpm/setup::remove-sdf *sim*
                                  (lambda (p)
                                    (cl-mpm/setup::plane-point-point-sdf
                                     p
                                     (cl-mpm/utils:vector-from-list (list ice-length datum 0d0))
                                     (cl-mpm/utils:vector-from-list (list (- ice-length cutback) offset 0d0))))
                                  :refine 3
                                  )
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
      (cl-mpm/setup::set-mass-filter *sim* density :proportion 1d-9)
      (setf (cl-mpm::sim-enable-fbar *sim*) fbar)
      (setf (cl-mpm/damage::sim-enable-length-localisation *sim*) t)
      (setf (cl-mpm::sim-allow-mp-split *sim*) t)
      (setf (cl-mpm::sim-max-split-depth *sim*) 6)
      ;; (setf (cl-mpm::sim-velocity-algorithm *sim*) :PIC)
      (setf (cl-mpm::sim-velocity-algorithm *sim*) :BLEND)



      (setf (cl-mpm:sim-dt *sim*)
            (* 0.5d0 (cl-mpm/setup:estimate-elastic-dt *sim*)))
      (setf (cl-mpm::sim-enable-damage *sim*) nil)
      (setf *run-sim* t)
      (defparameter *water-bc*
        (if t
             (cl-mpm/buoyancy::make-bc-buoyancy-clip
              *sim*
              datum
              water-density
              (lambda (pos datum) t)
              :visc-damping 1d-1
              )
             (cl-mpm/buoyancy::make-bc-buoyancy-body
              *sim*
              datum
              water-density
              (lambda (pos) t))))
      (cl-mpm:add-bcs-force-list *sim* *water-bc*)

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
           (* 1d9 0.01d0)
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
      (when (> offset 0d0) 
        (cl-mpm:add-bcs-force-list
         *sim*
         *floor-bc*
         ))
      ;; (cl-mpm:add-bcs-force-list
      ;;  *sim*
      ;;  *bc-erode*
      ;;  )
      (format t "MPs ~D~%" (length (cl-mpm:sim-mps *sim*)))
      )))


(defun run (&key (output-dir "./output/")
                 (ms 1d4)
                 )
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (uiop:ensure-all-directories-exist (list (uiop:merge-pathnames* output-dir)))
      (loop for f in (uiop:directory-files (uiop:merge-pathnames* output-dir)) do (uiop:delete-file-if-exists f))
      (with-open-file (stream (merge-pathnames output-dir "./timesteps.csv") :direction :output :if-exists :supersede)
        (format stream "steps,time,damage,plastic,energy,oobf,step-type,mass~%")))
    (cl-mpi:mpi-waitall)
    ;; (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
    ;;(cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~2,'0d_~5,'0d.vtk" rank 0)) *sim*)
    ;(cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_~2,'0d_~5,'0d.vtk" rank 0)) *sim*)
    ;(cl-mpm/output::save-vtk-nodes (merge-pathnames output-dir (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank 0)) *sim*)
    (defparameter *damping-water* 0d0)
    (let ((dt-scale 0.5d0)
          (visc-damping (cl-mpm/buoyancy::bc-viscous-damping *water-bc*))
          (fbar (cl-mpm::sim-enable-fbar *sim*))
          )
      (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 1d-2
               (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup:estimate-critical-damping *sim*)))
      ;;Fbar kinda bops our convergence
      (setf (cl-mpm::sim-enable-fbar *sim*) nil)
      ;;Extra viscous damping in the system is probably bad
      (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) 0d0)
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :oobf-crit 1d-2
       :energy-crit 1d-2
       :dt-scale dt-scale
       :substeps 50
       :conv-steps 1000
       :post-iter-step
       (lambda (i oobf energy)
         (cl-mpi:mpi-waitall)
         (cl-mpm/output:save-vtk (merge-pathnames output-dir (format nil "sim_conv_~2,'0d_~5,'0d.vtk" rank (+ 1 i))) *sim*)))
      (setf (cl-mpm::sim-damping-algorithm *sim*) :VISCOUS)
      (setf (cl-mpm:sim-damping-factor *sim*) 0d0)
      (setf (cl-mpm::sim-enable-fbar *sim*) fbar)
      (setf (cl-mpm/buoyancy::bc-viscous-damping *water-bc*) visc-damping)
      )
    (cl-mpm:iterate-over-mps 
      (cl-mpm:sim-mps *sim*) 
      (lambda (mp) 
        (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp))))

    (cl-mpm::iterate-over-mps
     (cl-mpm:sim-mps *sim*)
     (lambda (mp)
       (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))
    (setf (cl-mpm::sim-enable-damage *sim*) t)
    (setf (cl-mpm:sim-mass-scale *sim*) 1d0)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-6
             (cl-mpm/setup:estimate-critical-damping *sim*)))
    ;; (setf (cl-mpm/buoyancy::bc-enable *bc-erode*) t)
    (let* ((dt-scale 0.5)
           (substeps 0d0)
           (work 0d0)
           (oobf 0d0)
           (energy 0d0)
           (sim-state :accelerate)
           (accelerate-target-time 1d1)
           (accelerate-mass-scale ms)
           (collapse-target-time 1d0)
           (collapse-mass-scale 1d0)
           (criteria-energy 2d-2)
           (criteria-oobf 2d-2)
           (hist-factor 1.5)
           (target-time 2d0)
           (time 0d0)
           (crit-damp (cl-mpm/setup:estimate-critical-damping *sim*))
           )
      (when (= rank 0)
      (cl-mpm/output::save-simulation-parameters (merge-pathnames output-dir "settings.json")
                                     *sim*
                                     (list :dt target-time
                                           :criteria-energy criteria-energy
                                           :criteria-oobf criteria-oobf
                                           :ocean-height *water-height*
                                           :criteria-hist hist-factor 
                                           :delay-time *delay-time*
                                           )))
      (setf (cl-mpm::sim-mass-scale *sim*) accelerate-mass-scale
            target-time accelerate-target-time)
      (setf (cl-mpm:sim-dt *sim*) (* dt-scale (cl-mpm/setup:estimate-elastic-dt *sim*)))
      (setf substeps (ceiling target-time (cl-mpm:sim-dt *sim*)))
      (format t "Substeps ~D~%" substeps)
      (let ((damage-est 0d0)
            (mass-est 0d0))
        (loop for step from 0 below 1000
              while *run-sim*
              do
                 (setf damage-est
                       (cl-mpm/mpi:mpi-sum
                        (lparallel:pmap-reduce (lambda (mp)
                                                 (*
                                                  (cl-mpm/particle::mp-mass mp)
                                                  (cl-mpm/particle::mp-damage mp)))
                                               #'+ (cl-mpm:sim-mps *sim*)
                                               :initial-value 0d0))
                       )
                (setf mass-est 
                       (cl-mpm/mpi:mpi-sum
                        (lparallel:pmap-reduce (lambda (mp)
                                                 (cl-mpm/particle::mp-mass mp))
                                               #'+ (cl-mpm:sim-mps *sim*)
                                               :initial-value 0d0)))
                (let ((mp-count (cl-mpm/mpi:mpi-sum (float (length (cl-mpm:sim-mps *sim*)) 0d0))))
                  (when (= rank 0)
                    (format t "Step ~D~%" step)
                    (format t "Total MPs ~D~%" (round mp-count))
                    (with-open-file (stream (merge-pathnames output-dir "timesteps.csv") :direction :output :if-exists :append)
                      (format stream "~D,~f,~f,~f,~f,~f,~A,~f~%"
                              step
                              time
                              ;(cl-mpm::sim-time *sim*)
                              damage-est
                              0d0 ;; *data-plastic*
                              energy
                              oobf
                              sim-state
                              mass-est
                              ;; *data-mass*
                              ))))
                 ;; (setf work 0d0)
                 (cl-mpm/mpi::load-balance-algo *sim* :step-size 1d-2)
                 (setf energy 0d0
                       oobf 0d0
                       )
                 (time
                  (dotimes (i substeps)
                    (cl-mpm:update-sim *sim*)
                    (incf time (cl-mpm:sim-dt *sim*))
                    (incf oobf (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                    (incf energy (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                    (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))))
                 (setf
                  energy (/ energy substeps)
                  oobf (/ oobf substeps))

                 (if (= work 0d0)
                     (setf energy 0d0)
                     (setf energy (abs (/ energy work))))
                (when (or (>= energy (* criteria-energy hist-factor))
                          (>= oobf   (* criteria-oobf   hist-factor)))
                    (when (not (eq sim-state :collapse))
                       (setf sim-state :collapse)
                       (when (= rank 0)
                         (format t "Changed to collapse~%"))
                       (setf work 0d0)
                       )
                 )
                (when (and (< energy (/ criteria-energy hist-factor))
                           (< oobf   (/ criteria-oobf   hist-factor))) 
                  (progn
                       (when (not (eq sim-state :accelerate))
                         (when (= rank 0)
                           (format t "Changed to accelerate~%"))
                         (setf work 0d0)
                         (setf sim-state :accelerate)
                         (cl-mpm:iterate-over-mps
                          (cl-mpm:sim-mps *sim*)
                          (lambda (mp)
                            (cl-mpm/fastmaths::fast-zero (cl-mpm/particle:mp-velocity mp)))))))
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
                 (setf (cl-mpm:sim-damping-factor *sim*)
                       (* 1d-6 (sqrt (cl-mpm:sim-mass-scale *sim*)) crit-damp))
                 (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                   (when (= rank 0)
                     (format t "CFL dt estimate: ~f~%" dt-e)
                     (format t "CFL step count estimate: ~D~%" substeps-e))
                   (setf substeps substeps-e))
                 (swank.live:update-swank))))))

(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let* ((rank (cl-mpi:mpi-comm-rank))
        (cutback 000d0)
        (height (if (uiop:getenv "HEIGHT") (parse-float:parse-float (uiop:getenv "HEIGHT")) 400d0))
        (refine (if (uiop:getenv "REFINE") (parse-float:parse-float (uiop:getenv "REFINE")) 0.5))
        (fbar (if (uiop:getenv "FBAR") (string= (uiop:getenv "FBAR") "TRUE") nil))
        (ms (if (uiop:getenv "MS") (parse-float:parse-float (uiop:getenv "MS")) 1d4))
        (output-dir (merge-pathnames (format nil "/nobackup/rmvn14/paper-2/ice-buoyant-trial/output-~F_~F_~F_~D_~A/" height cutback refine (cl-mpi:mpi-comm-size) fbar)))
        )
    (format t "Outputting to - ~A~%" output-dir)
    (setup :refine refine
           :height height
           :mps 2
           :cutback cutback
           :fbar fbar)
    (when (typep *sim* 'cl-mpm/mpi::mpm-sim-mpi)
      (let* ((height (if (> (cl-mpi:mpi-comm-size) 1) 2 1))
             (dsize (ceiling (cl-mpi:mpi-comm-size) height)))
        (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize height 1)))
       (let ((mp (aref (cl-mpm:sim-mps *sim*) 0))) 
         (when (slot-exists-p mp 'cl-mpm/particle::local-length) 
           (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0))))) 
             (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))))
      (cl-mpm/mpi::setup-domain-bounds 
        *sim* 
        :domain-scaler (cl-mpm/mpi::origin-domain-scaler *sim*
                                                         *ice-length* 
                                                         *ice-height* 
                                                         nil))
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
    (run :output-dir output-dir)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  )

;(defparameter *output-directory* (merge-pathnames "/nobackup/rmvn14/paper-2/ice-buoyant/"))
;(defparameter *output-directory* (merge-pathnames "./output/"))
(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "16"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
;; (setup)
;; (format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;; (run :output-dir *output-directory*)



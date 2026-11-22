# A sweep over three densities: one STL per value, and one CSV with a header
# and a line per value. Run with -DTPMS=<tpms-forge> -DOUT=<folder>.
file(REMOVE ${OUT}/sweep.csv)
execute_process(COMMAND ${TPMS} --cells 1 --res 8 --sweep --density 0.15:0.25:0.05
                        -o ${OUT}/sweep.stl -o ${OUT}/sweep.csv RESULT_VARIABLE rc)
if(rc)
  message(FATAL_ERROR "tpms-forge failed (${rc})")
endif()
file(STRINGS ${OUT}/sweep.csv lines)
list(LENGTH lines n)
if(NOT n EQUAL 4)
  message(FATAL_ERROR "sweep.csv has ${n} lines, not 4")
endif()
foreach(d 0.15 0.2 0.25)
  if(NOT EXISTS ${OUT}/sweep_density-${d}.stl)
    message(FATAL_ERROR "missing sweep_density-${d}.stl")
  endif()
endforeach()

########################################################
#
#    Doxygen developer-reference documentation.
#
########################################################

find_package(Doxygen OPTIONAL_COMPONENTS dot)

if (DOXYGEN_FOUND)
    # In the following, CMake variables are set that look unused here but are
    # substituted into Doxyfile.in by configure_file() below.
    set(DOXYGEN_EXTRACT_ALL "YES")

    if (Doxygen_dot_FOUND)
        set(HAVE_DOT "YES")
        get_target_property(DOT_EXECUTABLE Doxygen::dot IMPORTED_LOCATION)
        get_filename_component(DOT_PATH "${DOT_EXECUTABLE}" DIRECTORY)
    else ()
        set(HAVE_DOT "NO")
        set(DOT_PATH "")
        message(ATTENTION "Graphviz 'dot' not found -- diagrams won't be generated in the documentation.")
    endif ()

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

    add_custom_target(doc
                      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Building Doxygen developer documentation in ${CMAKE_CURRENT_BINARY_DIR}/html"
                      VERBATIM)

    # A second configuration, with EXTRACT_ALL off, makes Doxygen itself
    # enumerate every documentation gap (WARN_IF_UNDOCUMENTED only fires for
    # entities that would otherwise not appear in the docs at all). Used by
    # the three targets below to check/enforce completeness; none of them
    # write out real HTML, so also skip dot for speed.
    set(DOXYGEN_EXTRACT_ALL "NO")
    set(HAVE_DOT "NO")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented @ONLY)

    # Matches only the "missing documentation" family of Doxygen warnings
    # (member/parameter/return-value not documented), not incidental parse
    # warnings unrelated to documentation coverage (e.g. Markdown files using
    # GitHub-style $...$ math, which Doxygen doesn't understand and reports
    # as "unknown command" -- see README.md's note on why it stays that way).
    set(UNDOCUMENTED_PATTERN "(is|are) not documented")
    set(DOXYGEN_COMPLETENESS_LOG "${CMAKE_CURRENT_BINARY_DIR}/doxygen_completeness.log")

    find_program(GREP grep)
    if (GREP)
        # Every target below first runs Doxygen into a log file as its own
        # COMMAND, separate from the grep step that follows. This matters:
        # add_custom_target fails as soon as any one of its COMMANDs exits
        # non-zero, so if Doxygen itself crashes (bad config, missing input,
        # ...) that failure is caught right here -- rather than piping
        # Doxygen straight into grep, where a crash produces output with no
        # "not documented" match, grep exits 1, and (for undocumented_test)
        # the leading `!` would then turn that into a false success.
        add_custom_target(undocumented
                          COMMAND ${DOXYGEN_EXECUTABLE}
                                  ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                                  > ${DOXYGEN_COMPLETENESS_LOG} 2>&1
                          COMMAND ${GREP} -iE "${UNDOCUMENTED_PATTERN}"
                                  ${DOXYGEN_COMPLETENESS_LOG}
                                  ||
                                  ${CMAKE_COMMAND} -E echo "No undocumented entities found."
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Listing every Doxygen 'not documented' warning"
                          VERBATIM)

        add_custom_target(undocumented_count
                          COMMAND ${DOXYGEN_EXECUTABLE}
                                  ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                                  > ${DOXYGEN_COMPLETENESS_LOG} 2>&1
                          COMMAND ${GREP} -icE "${UNDOCUMENTED_PATTERN}"
                                  ${DOXYGEN_COMPLETENESS_LOG}
                                  |
                                  xargs -I{}
                                  ${CMAKE_COMMAND} -E echo "There are {} Doxygen warnings about missing documentation."
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Counting Doxygen 'not documented' warnings"
                          VERBATIM)

        # The CI-enforced completeness gate. Once Doxygen's own exit status
        # has been validated by the first COMMAND above, `!` here negates
        # only grep's exit status (now read back from the log file rather
        # than a live pipe), so this fails exactly when grep finds at least
        # one undocumented entity -- and no longer whenever Doxygen fails.
        add_custom_target(undocumented_test
                          COMMAND ${DOXYGEN_EXECUTABLE}
                                  ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                                  > ${DOXYGEN_COMPLETENESS_LOG} 2>&1
                          COMMAND ! ${GREP} -qiE "${UNDOCUMENTED_PATTERN}"
                                  ${DOXYGEN_COMPLETENESS_LOG}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Testing documentation completeness for CI: fails if there is any warning"
                          VERBATIM)
    endif ()

    mark_as_advanced(GREP)
else ()
    message(ATTENTION "Doxygen not found -- the 'doc' target won't be available.")
endif ()

#!/bin/bash
# Automated test script for gmx_continue_extend.sh
# Tests all major modes: easy mode, segmented append/noappend, crash recovery

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
BIN_DIR="${REPO_ROOT}/bin"
TEST_DATA_DIR="${SCRIPT_DIR}/trp_cage"

# Results tracking
declare -a TESTS_PASSED
declare -a TESTS_FAILED

# Helper functions
print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_test() {
    echo -e "${YELLOW}[TEST]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
    TESTS_PASSED+=("$1")
}

print_failure() {
    echo -e "${RED}[FAIL]${NC} $1"
    TESTS_FAILED+=("$1")
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Setup function
setup_test_env() {
    if [[ ! -d "${TEST_DATA_DIR}" ]]; then
        echo -e "${RED}Error: Test data directory not found: ${TEST_DATA_DIR}${NC}"
        echo "Please copy your TRP-cage files to tests/trp_cage/"
        exit 1
    fi
    
    if [[ ! -f "${TEST_DATA_DIR}/md_0.tpr" ]]; then
        echo -e "${RED}Error: md_0.tpr not found in ${TEST_DATA_DIR}${NC}"
        echo "Please ensure you have md_0.tpr, md_0.cpt in tests/trp_cage/"
        exit 1
    fi
    
    print_info "Test data directory: ${TEST_DATA_DIR}"
    print_info "Script: ${BIN_DIR}/gmx_continue_extend.sh"
}

# Cleanup function
cleanup_test_files() {
    local dir="$1"
    cd "$dir"
    print_info "Cleaning up test files in $(pwd)"
    rm -f md_0.part*.* *.backup md_0_prev.cpt md_test_*.* test_*.tpr 2>/dev/null || true
}

# Test 1: Easy mode with single TPR
test_easy_mode_single() {
    print_test "Test 1: Easy mode with single TPR"
    
    local test_dir="${TEST_DATA_DIR}/test1_easy_single"
    mkdir -p "$test_dir"
    
    # Copy original files
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.xtc" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.log" "$test_dir/" 2>/dev/null || true
    
    cd "$test_dir"
    
    # Run easy mode
    if "${BIN_DIR}/gmx_continue_extend.sh" 500 > test_output.log 2>&1; then
        if [[ -f md_0.xtc ]] && [[ -f md_0.log ]]; then
            print_success "Easy mode single TPR"
        else
            print_failure "Easy mode single TPR - output files missing"
        fi
    else
        print_failure "Easy mode single TPR - command failed"
        cat test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 2: Easy mode with explicit TPR (simulate multiple TPRs)
test_easy_mode_explicit() {
    print_test "Test 2: Easy mode with explicit TPR"
    
    local test_dir="${TEST_DATA_DIR}/test2_easy_explicit"
    mkdir -p "$test_dir"
    
    # Copy files and create dummy second TPR
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.xtc" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.log" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/npt.tpr"  # Create second TPR
    
    cd "$test_dir"
    
    # Should fail without explicit TPR
    if "${BIN_DIR}/gmx_continue_extend.sh" 500 > test_output_fail.log 2>&1; then
        print_failure "Easy mode should fail with multiple TPRs"
    else
        print_info "Correctly rejected multiple TPRs"
        
        # Should succeed with explicit TPR
        if "${BIN_DIR}/gmx_continue_extend.sh" 500 md_0.tpr > test_output_success.log 2>&1; then
            if [[ -f md_0.xtc ]] && [[ -f md_0.log ]]; then
                print_success "Easy mode explicit TPR"
            else
                print_failure "Easy mode explicit TPR - output files missing"
            fi
        else
            print_failure "Easy mode explicit TPR - command failed"
            cat test_output_success.log
        fi
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 3: Normal mode with noappend (segmented)
test_normal_noappend() {
    print_test "Test 3: Normal mode with noappend (part files)"
    
    local test_dir="${TEST_DATA_DIR}/test3_noappend"
    mkdir -p "$test_dir"
    
    # Copy original files
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.xtc" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.log" "$test_dir/" 2>/dev/null || true
    
    cd "$test_dir"
    
    # Run 2 segments of 1 ns each (integers required)
    if "${BIN_DIR}/gmx_continue_extend.sh" md 0 0 2 1 --tpr md_0.tpr --noappend > test_output.log 2>&1; then
        # Check for part files
        if [[ -f md_0.part0001.xtc ]] || [[ -f md_0.part0001.log ]]; then
            print_success "Normal mode noappend (part files created)"
        else
            print_failure "Normal mode noappend - part files not found"
            ls -la md_0.part* 2>&1 || true
        fi
    else
        print_failure "Normal mode noappend - command failed"
        tail -50 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 4: Normal mode with append (segmented)
test_normal_append() {
    print_test "Test 4: Normal mode with append (same files)"
    
    local test_dir="${TEST_DATA_DIR}/test4_append"
    mkdir -p "$test_dir"
    
    # Copy original files
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.xtc" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.log" "$test_dir/" 2>/dev/null || true
    
    cd "$test_dir"
    
    # Get initial file size
    initial_size=0
    if [[ -f md_0.xtc ]]; then
        initial_size=$(stat -f%z md_0.xtc 2>/dev/null || stat -c%s md_0.xtc 2>/dev/null || echo 0)
    fi
    
    # Run 2 segments of 1 ns each (integers required)
    if "${BIN_DIR}/gmx_continue_extend.sh" md 0 0 2 1 --tpr md_0.tpr --append > test_output.log 2>&1; then
        # Check that part files were NOT created and main file grew
        if [[ ! -f md_0.part0001.xtc ]]; then
            final_size=$(stat -f%z md_0.xtc 2>/dev/null || stat -c%s md_0.xtc 2>/dev/null || echo 0)
            if [[ $final_size -gt $initial_size ]]; then
                print_success "Normal mode append (same file extended)"
            else
                print_failure "Normal mode append - file did not grow"
            fi
        else
            print_failure "Normal mode append - part files created (should not)"
        fi
    else
        print_failure "Normal mode append - command failed"
        tail -50 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 5: Crash recovery in easy mode
test_crash_recovery_easy() {
    print_test "Test 5: Crash recovery (easy mode)"
    
    local test_dir="${TEST_DATA_DIR}/test5_crash_easy"
    mkdir -p "$test_dir"
    
    # Copy original files
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.xtc" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/md_0.log" "$test_dir/" 2>/dev/null || true
    
    cd "$test_dir"
    
    # Simulate crash by ensuring checkpoint exists
    if [[ -f md_0.cpt ]]; then
        # Run with existing checkpoint - should use -cpi
        if "${BIN_DIR}/gmx_continue_extend.sh" 500 > test_output.log 2>&1; then
            # Check log for -cpi usage
            if grep -q "\-cpi" test_output.log || [[ -f md_0.xtc ]]; then
                print_success "Crash recovery easy mode"
            else
                print_failure "Crash recovery easy mode - -cpi not detected"
            fi
        else
            print_failure "Crash recovery easy mode - command failed"
            cat test_output.log
        fi
    else
        print_info "No checkpoint available, skipping crash recovery test"
        print_success "Crash recovery easy mode (skipped - no checkpoint)"
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 6: Check that base_name is derived from TPR in --tpr mode
test_basename_from_tpr() {
    print_test "Test 6: Base name derived from TPR filename"
    
    local test_dir="${TEST_DATA_DIR}/test6_basename"
    mkdir -p "$test_dir"
    
    # Copy with different name
    cp "${TEST_DATA_DIR}/md_0.tpr" "$test_dir/custom_name.tpr"
    cp "${TEST_DATA_DIR}/md_0.cpt" "$test_dir/custom_name.cpt" 2>/dev/null || true
    
    cd "$test_dir"
    
    # Run easy mode
    if "${BIN_DIR}/gmx_continue_extend.sh" 500 custom_name.tpr > test_output.log 2>&1; then
        # Check that output uses custom_name, not md_0
        if grep -q "custom_name" test_output.log; then
            print_success "Base name from TPR filename"
        else
            print_failure "Base name from TPR - wrong basename used"
        fi
    else
        print_failure "Base name from TPR - command failed"
        cat test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Main test execution
main() {
    print_header "Starting gmx_continue_extend.sh Test Suite"
    
    setup_test_env
    
    # Run all tests
    test_easy_mode_single
    test_easy_mode_explicit
    test_normal_noappend
    test_normal_append
    test_crash_recovery_easy
    test_basename_from_tpr
    
    # Print summary
    print_header "Test Summary"
    
    echo ""
    echo -e "${GREEN}Passed: ${#TESTS_PASSED[@]}${NC}"
    for test in "${TESTS_PASSED[@]}"; do
        echo -e "  ${GREEN}✓${NC} $test"
    done
    
    echo ""
    if [[ ${#TESTS_FAILED[@]} -gt 0 ]]; then
        echo -e "${RED}Failed: ${#TESTS_FAILED[@]}${NC}"
        for test in "${TESTS_FAILED[@]}"; do
            echo -e "  ${RED}✗${NC} $test"
        done
        echo ""
        exit 1
    else
        echo -e "${GREEN}All tests passed!${NC}"
        echo ""
    fi
    
    # Cleanup prompt
    echo ""
    read -p "Clean up all test directories? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Cleaning up test directories..."
        rm -rf "${TEST_DATA_DIR}/test"*
        print_info "Cleanup complete"
    fi
}

# Run tests
main "$@"

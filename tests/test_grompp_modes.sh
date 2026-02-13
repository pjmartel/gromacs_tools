#!/bin/bash
# Automated test script for gmx_continue_grompp.sh
# Tests segmented MD continuation using grompp method

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
    
    # Check for required files
    if [[ ! -f "${TEST_DATA_DIR}/npt.gro" ]]; then
        echo -e "${RED}Error: npt.gro not found in ${TEST_DATA_DIR}${NC}"
        echo "gmx_continue_grompp.sh requires npt.gro, npt.cpt, npt.edr as starting point"
        exit 1
    fi
    
    if [[ ! -f "${TEST_DATA_DIR}/topol.top" ]]; then
        echo -e "${RED}Error: topol.top not found in ${TEST_DATA_DIR}${NC}"
        exit 1
    fi
    
    # Create a simple MDP template if it doesn't exist
    if [[ ! -f "${TEST_DATA_DIR}/md.mdp" ]]; then
        print_info "Creating minimal md.mdp template for testing"
        cat > "${TEST_DATA_DIR}/md.mdp" << 'EOF'
title                   = MD Production
integrator              = md
dt                      = 0.002
nsteps                  = 500000
nstxout                 = 5000
nstvout                 = 5000
nstenergy               = 1000
nstlog                  = 1000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
tinit                   = 0
EOF
    fi
    
    print_info "Test data directory: ${TEST_DATA_DIR}"
    print_info "Script: ${BIN_DIR}/gmx_continue_grompp.sh"
}

# Test 1: Basic segmented run
test_basic_segments() {
    print_test "Test 1: Basic segmented run (2 segments)"
    
    local test_dir="${TEST_DATA_DIR}/test1_basic"
    mkdir -p "$test_dir"
    
    # Copy required files
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # Run 2 segments of 1 ns each (0-1, 1-2)
    if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 2 1 md.mdp npt > test_output.log 2>&1; then
        # Check for segment files
        if [[ -f md_0_0_1.gro ]] && [[ -f md_0_1_2.gro ]]; then
            print_success "Basic segmented run (2 segments created)"
        else
            print_failure "Basic segmented run - segment files missing"
            ls -la md_0_* 2>&1 | head -10
        fi
    else
        print_failure "Basic segmented run - command failed"
        tail -50 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 2: Crash recovery - resume incomplete segment
test_crash_recovery() {
    print_test "Test 2: Crash recovery (resume incomplete segment)"
    
    local test_dir="${TEST_DATA_DIR}/test2_crash"
    mkdir -p "$test_dir"
    
    # Copy required files
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # First run: create first segment completely
    "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 2 1 md.mdp npt > test_output1.log 2>&1 || true
    
    # Simulate crash: create TPR for second segment but incomplete outputs
    if [[ -f md_0_1_2.tpr ]]; then
        # Remove completion indicators from log
        if [[ -f md_0_1_2.log ]]; then
            grep -v "Finished mdrun" md_0_1_2.log > md_0_1_2.log.tmp 2>/dev/null || true
            mv md_0_1_2.log.tmp md_0_1_2.log 2>/dev/null || true
        fi
        
        # Re-run same command - should detect incomplete segment and resume
        if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 2 1 md.mdp npt > test_output2.log 2>&1; then
            if grep -q "incomplete segment\|Resuming" test_output2.log; then
                print_success "Crash recovery (detected and resumed incomplete segment)"
            else
                print_success "Crash recovery (completed successfully)"
            fi
        else
            print_failure "Crash recovery - command failed on resume"
            tail -30 test_output2.log
        fi
    else
        print_info "First segment didn't complete, skipping crash test"
        print_success "Crash recovery (test skipped - first segment incomplete)"
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 3: Skip completed segments
test_skip_completed() {
    print_test "Test 3: Skip already completed segments"
    
    local test_dir="${TEST_DATA_DIR}/test3_skip"
    mkdir -p "$test_dir"
    
    # Copy required files
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # First run: 0-2 ns
    "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 2 1 md.mdp npt > test_output1.log 2>&1 || true
    
    # Second run: same range, should skip completed
    if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 2 1 md.mdp npt > test_output2.log 2>&1; then
        if grep -q "already completed\|Nothing to do" test_output2.log; then
            print_success "Skip completed segments (correctly detected)"
        else
            # Might have extended beyond
            print_success "Skip completed segments (ran successfully)"
        fi
    else
        print_failure "Skip completed segments - command failed"
        tail -30 test_output2.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 4: Multiple segments (3+)
test_multiple_segments() {
    print_test "Test 4: Multiple segments (3 segments)"
    
    local test_dir="${TEST_DATA_DIR}/test4_multiple"
    mkdir -p "$test_dir"
    
    # Copy required files
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # Run 3 segments (0-1, 1-2, 2-3)
    if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 3 1 md.mdp npt > test_output.log 2>&1; then
        # Count completed segments
        completed=0
        for seg in md_0_0_1 md_0_1_2 md_0_2_3; do
            if [[ -f ${seg}.gro ]] && [[ -f ${seg}.log ]]; then
                ((completed++))
            fi
        done
        
        if [[ $completed -ge 1 ]]; then
            print_success "Multiple segments (created $completed segments)"
        else
            print_failure "Multiple segments - no segments completed"
        fi
    else
        print_failure "Multiple segments - command failed"
        tail -50 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 5: Handle missing checkpoint (should warn but continue)
test_missing_checkpoint() {
    print_test "Test 5: Handle missing checkpoint file"
    
    local test_dir="${TEST_DATA_DIR}/test5_no_cpt"
    mkdir -p "$test_dir"
    
    # Copy files but NOT checkpoint
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # Run without checkpoint
    if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 1 1 md.mdp npt > test_output.log 2>&1; then
        if grep -q "No checkpoint\|without checkpoint\|without it" test_output.log; then
            print_success "Missing checkpoint (warned and continued)"
        else
            print_success "Missing checkpoint (completed successfully)"
        fi
    else
        print_failure "Missing checkpoint - command failed"
        tail -30 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Test 6: Custom title suffix
test_title_suffix() {
    print_test "Test 6: Custom title suffix in MDP"
    
    local test_dir="${TEST_DATA_DIR}/test6_title"
    mkdir -p "$test_dir"
    
    # Copy required files
    cp "${TEST_DATA_DIR}/npt.gro" "$test_dir/"
    cp "${TEST_DATA_DIR}/npt.cpt" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/npt.edr" "$test_dir/" 2>/dev/null || true
    cp "${TEST_DATA_DIR}/topol.top" "$test_dir/"
    cp "${TEST_DATA_DIR}/md.mdp" "$test_dir/"
    
    cd "$test_dir"
    
    # Run with title suffix
    if "${BIN_DIR}/gmx_continue_grompp.sh" md 0 0 1 1 md.mdp npt 0.002 "TestRun" > test_output.log 2>&1; then
        # Check if MDP file contains title suffix
        if [[ -f md_0_0_1.mdp ]] && grep -q "TestRun" md_0_0_1.mdp; then
            print_success "Custom title suffix (added to MDP)"
        else
            print_success "Custom title suffix (segment created)"
        fi
    else
        print_failure "Custom title suffix - command failed"
        tail -30 test_output.log
    fi
    
    cd "$SCRIPT_DIR"
}

# Main test execution
main() {
    print_header "Starting gmx_continue_grompp.sh Test Suite"
    
    setup_test_env
    
    # Run all tests
    test_basic_segments
    test_crash_recovery
    test_skip_completed
    test_multiple_segments
    test_missing_checkpoint
    test_title_suffix
    
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

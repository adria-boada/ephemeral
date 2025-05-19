#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# predefined_models.py
#
# 02 de març 2025  <adria@molevol-OptiPlex-9020>

from . import Models_2D

# AVAILABLE TESTING/FITTING/OPTIMIZING MODELS
# -------------------------------------------
models = {
    #> Standard neutral model, populations never diverge.
    #>"no_divergence": {
    #>    "func": Models_2D.no_divergence,
    #>    "param": 1,
    #>    "labels": "nu",
    #>},
    # Split into two populations, no migration.
    "no_mig": {
        "func": Models_2D.no_mig,
        "param": 3,
        "labels": "nu1, nu2, T",
        "lower": [1e-3, 1e-3, 0],
        "upper": [100,   100, 40],
    },
    # Split into two populations, with symmetric migration.
    "sym_mig": {
        "func": Models_2D.sym_mig,
        "param": 4,
        "labels": "nu1, nu2, m, T",
        "lower": [1e-3, 1e-3,   0, 0],
        "upper": [100,   100, 1e4, 40],
    },
    # Split into two populations, with different migration rates.
    "asym_mig": {
        "func": Models_2D.asym_mig,
        "param": 5,
        "labels": "nu1, nu2, m12, m21, T",
        "lower": [1e-3, 1e-3,   0,   0, 0],
        "upper": [100,   100, 1e4, 1e4, 40],
    },
    # Split with symmetric migration followed by isolation.
    "anc_sym_mig": {
        "func": Models_2D.anc_sym_mig,
        "param": 5,
        "labels": "nu1, nu2, m, T1, T2",
        "lower": [1e-3, 1e-3,   0, 0, 0],
        "upper": [100,   100, 1e4, 40, 40],
    },
    # Split with asymmetric migration followed by isolation.
    "anc_asym_mig": {
        "func": Models_2D.anc_asym_mig,
        "param": 6,
        "labels": "nu1, nu2, m12, m21, T1, T2",
        "lower": [1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [100,   100, 1e4, 1e4, 40, 40],
    },
    # Split with no gene flow, followed by a
    # period of symmetrical gene flow.
    "sec_contact_sym_mig": {
        "func": Models_2D.sec_contact_sym_mig,
        "param": 5,
        "labels": "nu1, nu2, m, T1, T2",
        "lower": [1e-3, 1e-3,   0, 0, 0],
        "upper": [100,   100, 1e4, 40, 40],
    },
    # Split with no gene flow, followed by a
    # period of asymmetrical gene flow.
    "sec_contact_asym_mig": {
        "func": Models_2D.sec_contact_asym_mig,
        "param": 6,
        "labels": "nu1, nu2, m12, m21, T1, T2",
        "lower": [1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [100,   100, 1e4, 1e4, 40, 40],
    },
    
    #-MODELS INVOLVING SIZE CHANGES-#

    # Split with no migration, then size change
    # with no migration.
    "no_mig_size": {
        "func": Models_2D.no_mig_size,
        "param": 6,
        "labels": "nu1a, nu2a, nu1b, nu2b, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3, 0, 0],
        "upper": [ 100,  100,  100,  100, 40, 40],
    },
    # Split with symmetric migration, then size change
    # with symmetric migration.
    "sym_mig_size": {
        "func": Models_2D.sym_mig_size,
        "param": 7,
        "labels": "nu1a, nu2a, nu1b, nu2b, m, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 40, 40],
    },
    # Split with different migration rates, then
    # size change with different migration rates.
    "asym_mig_size": {
        "func": Models_2D.asym_mig_size,
        "param": 8,
        "labels": "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 1e4, 40, 40],
    },
    # Split with symmetrical gene flow, followed by size change
    # with no gene flow.  
    "anc_sym_mig_size": {
        "func": Models_2D.anc_sym_mig_size,
        "param": 7,
        "labels": "nu1a, nu2a, nu1b, nu2b, m, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 40, 40],
    },
    # Split with asymmetrical gene flow, followed by size change
    # with no gene flow.
    "anc_asym_mig_size": {
        "func": Models_2D.anc_asym_mig_size,
        "param": 8,
        "labels": "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 1e4, 40, 40],
    },
    # Split with no gene flow, followed by size change
    # with symmetrical gene flow.
    "sec_contact_sym_mig_size": {
        "func": Models_2D.sec_contact_sym_mig_size,
        "param": 7,
        "labels": "nu1a, nu2a, nu1b, nu2b, m, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 40, 40],
    },
    # Split with no gene flow, followed by size change
    # with asymmetrical gene flow.
    "sec_contact_asym_mig_size": {
        "func": Models_2D.sec_contact_asym_mig_size,
        "param": 8,
        "labels": "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 1e4, 40, 40],
    },

    #-TWO EPOCH SPLIT WITH CHANGING MIG. RATE-#

    # Split into two populations, with symmetric migration.
    # A second period of symmetric migration occurs,
    # but can be a different rate. Pop size is same.
    "sym_mig_twoepoch": {
        "func": Models_2D.sym_mig_twoepoch,
        "param": 6,
        "labels": "nu1, nu2, m1, m2, T1, T2",
        "lower": [1e-3, 1e-3,   0,   0, 0, 0],
        "upper": [ 100,  100, 1e4, 1e4, 40, 40],
    },
    # Split into two populations, with different migration rates.
    # A second period of asymmetric migration occurs,
    # but can be at different rates. Pop size is same.
    "asym_mig_twoepoch": {
        "func": Models_2D.asym_mig_twoepoch,
        "param": 8,
        "labels": "nu1, nu2, m12a, m21a, m12b, m21b, T1, T2",
        "lower": [1e-3, 1e-3,   0,   0,   0,   0, 0, 0],
        "upper": [ 100,  100, 1e4, 1e4, 1e4, 1e4, 40, 40],
    },

    #-THREE EPOCH WITH SECONDARY CONTACT AND ISOLATION-#

    # Split with no gene flow, followed by a period of
    # symmetrical gene flow, then isolation.
    "sec_contact_sym_mig_three_epoch": {
        "func": Models_2D.sec_contact_sym_mig_three_epoch,
        "param": 6,
        "labels": "nu1, nu2, m, T1, T2, T3",
        "lower": [1e-3, 1e-3,   0, 0, 0, 0],
        "upper": [ 100,  100, 1e4, 40, 40, 40],
    },
    # Split with no gene flow, followed by a period of
    # asymmetrical gene flow, then isolation.
    "sec_contact_asym_mig_three_epoch": {
        "func": Models_2D.sec_contact_asym_mig_three_epoch,
        "param": 7,
        "labels": "nu1, nu2, m12, m21, T1, T2, T3",
        "lower": [1e-3, 1e-3,   0,   0, 0, 0, 0],
        "upper": [ 100,  100, 1e4, 1e4, 40, 40, 40],
    },
    # Split with no gene flow, followed by size change with
    # symmetrical gene flow, then isolation.
    "sec_contact_sym_mig_size_three_epoch": {
        "func": Models_2D.sec_contact_sym_mig_size_three_epoch,
        "param": 8,
        "labels": "nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0, 0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 40, 40, 40],
    },
    # Split with no gene flow, followed by size change with
    # asymmetrical gene flow, then isolation.
    "sec_contact_asym_mig_size_three_epoch": {
        "func": Models_2D.sec_contact_asym_mig_size_three_epoch,
        "param": 9,
        "labels": "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3",
        "lower": [1e-3, 1e-3, 1e-3, 1e-3,   0,   0, 0, 0, 0],
        "upper": [ 100,  100,  100,  100, 1e4, 1e4, 40, 40, 40],
    },
}


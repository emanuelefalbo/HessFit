from xtb_screening import run_xtb

result = run_xtb(
    "conformers/conf_000.xyz",
    "xtb_screening/conf_000"
)

print("Success:", result.success)
print("Energy:", result.energy)
print("Optimized:", result.optimized_xyz)
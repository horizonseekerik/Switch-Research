with open("/mnt/c/Users/hp/.gemini/antigravity/brain/5f75e54d-9002-4d83-af37-5cbe8ca658e8/.system_generated/tasks/task-5340.log") as f:
    for line in f:
        line_str = line.strip()
        if any(line_str.startswith(p) for p in ["4.20", "4.28", "4.35", "4.42", "L_c"]):
            print(line_str)

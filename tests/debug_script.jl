

macro kwfunc_from_dict(fname_expr, dict_expr, body)
    # fname_expr is e.g. :(myfunc(x, y))
    # dict_expr is the constant dictionary expression, e.g. :D
    # body is the function body expression

    # Evaluate dict_expr at macro expansion time to get keys and values
    dict_val = eval(dict_expr)
    if !(dict_val isa AbstractDict)
        error("Second argument must be a constant dictionary")
    end

    # Extract positional args and function name
    fname = nothing
    pos_args = []

    if fname_expr.head == :call
        fname = fname_expr.args[1]
        pos_args = fname_expr.args[2:end]
    else
        fname = fname_expr
    end

    # Build keyword argument expressions with defaults from dict_val
    kw_args = []
    for k in keys(dict_val)
        # keys can be Symbol or String, convert to Symbol if needed
        sym = k isa Symbol ? k : Symbol(k)
        val = dict_val[k]
        push!(kw_args, :( $(sym) = $(val) ))
    end

    # Compose function signature: positional args + keyword args
    fn_signature = Expr(:call, fname, pos_args..., :(;), kw_args...)

    # Compose full function expression
    fn_expr = Expr(:function, fn_signature, body)

    return fn_expr
end

D  = Dict(:a=>2.45,:b=>"dfg",:c=>234.5)
names = collect(keys(D))
vals = collect(values(D))


@kwfunc_from_dict(fill_dict(a),:D,2+2)   

Dict(
    "title" => @htl("Computational Data Analysis <strong>class website</strong>"),

    # # add a disclaimer to the course webpage. Remove it if you dont want to include it.
    "disclaimer" => md"""
    This website is designed to provide a comprehensive introduction to computational data analysis using the Julia programming language and the Pluto.jl interactive environment. If you are in the class, or simply want to learn a thing or two, feel free to look around!
    """,

    # Highlights the key features of your class to make it more engaging. Remove it if you dont want to include it.
    "highlights" => [
        Dict("name" => "Fundamental Characteristics of Data", 
             "text" => md"This explains how to describe, clean, and organize a collection of numbers. Representing data in a meaningful way is the first step in any data analysis project.",
             "img" => "https://github.com/quentingageclark-a11y/plutoimages/blob/main/stephen-dawson-qwtCeJ5cLYs-unsplash.jpg?raw=true"
        ),
        Dict("name" => "Signal Analysis",
             "text" => md"""
             We have a collection of lots of data over time and space. With the help of signal analysis, we can work to clean this data and extract meaningful information and patterns from background noise.
             """,
             "img" => "https://github.com/quentingageclark-a11y/plutoimages/blob/main/aedrian-salazar-e1e2JiEAd70-unsplash.jpg?raw=true"
             ),
        Dict("name" => "Numerical Analysis- breaking down reality into data points; what are the gaps?",
             "text" => md"""
             Numerical analysis is the study of algorithms for the problems of continuous mathematics (as distinguished from discrete mathematics).
             It is a branch of mathematics that deals with the development and analysis of numerical methods for solving mathematical problems.
             """,
             "img" => "https://user-images.githubusercontent.com/6933510/136203632-29ce0a96-5a34-46ad-a996-de55b3bcd380.png"
        )
    ]
)
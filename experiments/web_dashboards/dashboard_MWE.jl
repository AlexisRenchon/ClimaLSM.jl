using JSServe
import JSServe.TailwindDashboard as D
using WGLMakie

n = Observable(1)
function button_dashboard(button)
    on(button) do click
        n[] += 1
    end
    return n
end

App() do
    button = D.Button("click")
    n = button_dashboard(button)
    return DOM.div(
                   D.Card(button), D.Card(n)
                  )
end









































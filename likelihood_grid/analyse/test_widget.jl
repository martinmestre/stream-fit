using Interact, DynamicalSystems, Plots
using Mux, WebIO
using Blink

function app(req)
    @manipulate for px=0:0.005:0.3, py=0:0.005:0.3, plane = Dict("x=0" => 1, "y=0" => 2)
        psos = poincaresos(
            Systems.henonheiles(), (plane, 0), 10000; u0 = [0., 0., px, py])

            vbox(
                scatter(psos[:, [2,1][plane]], psos[:, [4,3][plane]],
                    linewidth=0.1, markersize=2, markerstrokewidth=0,
                    xlimit=[-1/3, 1/3], ylimit=[-1/3, 1/3],
                    markercolor="#F35700"
                ))
    end
end

# With Mux-WebIO (remember to open local web url):
#WebIO.webio_serve(page("/", app), 8004)

# With Blink (Electron terminal):
w = window()
body!(w, app(8010))
